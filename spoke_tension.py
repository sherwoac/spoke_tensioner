#!/usr/bin/env python3
"""
Spoke tension meter + wheel equaliser.

    pip install sounddevice numpy matplotlib
    python spoke_tension.py --device 0            # live
    python spoke_tension.py --file pluck.wav      # analyse one recording

Keys (in the plot window):
    space   start / stop listening (starts paused)
    enter / right arrow   next spoke          left arrow   previous spoke
    1-9 then enter        jump to spoke N     (type digits, enter to commit)
    c       calibrate turn gain: press once, turn the current nipple +1/4 turn, press again
    l       calibrate length: press, pluck with crossing partner held OFF (full length),
            then pluck normally; the effective vibrating length is computed and applied
    u       undo last reading       s   save session now
    q / esc / close window / Ctrl-C   quit cleanly (session is saved)

Session is auto-saved to --session (default wheel_session.json) after every reading.
"""

import argparse
import json
import os
import signal
import sys
import time

import numpy as np
import sounddevice as sd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

G = 9.80665
STEEL_DENSITY = 7900.0
THREAD_PITCH_MM = 0.4536          # 14 G spoke thread, 56 TPI
TONE_WIN_S = 0.15                 # sliding window for the tone detector
STEADY_WINDOWS = 2                # consecutive agreeing windows = a pluck


# ------------------------------------------------------------------ physics
def mass_per_metre(diameter_mm):
    r = diameter_mm / 2000.0
    return STEEL_DENSITY * np.pi * r * r


def freq_to_kgf(f, length_m, mu):
    return 4.0 * length_m ** 2 * f ** 2 * mu / G


def kgf_to_freq(kgf, length_m, mu):
    return np.sqrt(kgf * G / (4.0 * length_m ** 2 * mu))


# ---------------------------------------------------------- pitch detection
def bandpass_fft(x, fs, flo, fhi):
    n = len(x)
    spec = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, 1.0 / fs)
    width = 0.08 * (fhi - flo)
    mask = np.clip((freqs - (flo - width)) / width, 0, 1) * \
           np.clip(((fhi + width) - freqs) / width, 0, 1)
    return np.fft.irfft(spec * mask, n)


SEGMENT_RATIO = [1.46]   # f_segment / f_full-length = L_full / L_segment; set by main()


def spectral_peak(x, fs, fmin, fmax):
    """Harmonic-weighted spectral estimate of the fundamental inside [fmin, fmax].
    Returns (f, prominence, (freqs, mags) up to 2 kHz)."""
    n = len(x)
    nfft = 1 << int(np.ceil(np.log2(n * 8)))
    spec = np.abs(np.fft.rfft(x * np.hanning(n), nfft))
    df = fs / nfft
    keep = int(2000 / df)
    view = (np.arange(keep) * df, spec[:keep])
    lo, hi = int(fmin / df), int(fmax / df)
    if hi - lo < 5:
        return None, 0.0, view
    band = spec[lo:hi]
    floor = np.median(band) + 1e-12
    pk = np.where((band[1:-1] >= band[:-2]) & (band[1:-1] > band[2:]) &
                  (band[1:-1] > 2.0 * floor))[0] + 1 + lo
    if len(pk) == 0:
        return None, 0.0, view

    def mag_at(f):
        i = int(round(f / df)); w = max(1, int(0.015 * f / df))
        a, b = max(0, i - w), min(len(spec), i + w + 1)
        return float(np.max(spec[a:b])) if b > a else 0.0

    scored = []
    for i in pk:
        f = i * df
        scored.append((spec[i] + 0.7 * mag_at(2 * f) + 0.4 * mag_at(3 * f), i))
    scored.sort(reverse=True)
    best = scored[0][1]
    # loose-spoke guard: if the runner-up is a strong peak at SEGMENT_RATIO x the winner,
    # the winner is the full-length mode and the runner-up is the real segment tone.
    r = SEGMENT_RATIO[0]
    if r > 1.15:                      # only meaningful when the segment is clearly shorter than the spoke
        for sc, j in scored[1:4]:
            if sc > 0.4 * scored[0][0] and abs(j / best - r) / r < 0.06:
                best = j
                break
    # octave-down guard: if there is a strong peak near HALF the winner's frequency,
    # the winner is a harmonic / upper mode and the lower peak is the fundamental.
    for _ in range(2):
        half = best / 2
        if half * df < fmin:
            break
        cands = [j for j in pk if abs(j - half) / half < 0.08 and spec[j] >= 0.5 * spec[best]]
        if not cands:
            break
        best = max(cands, key=lambda j: spec[j])
    i = best
    y0, y1, y2 = spec[i - 1], spec[i], spec[i + 1]
    d = y0 - 2 * y1 + y2
    off = float(np.clip(0.5 * (y0 - y2) / d if abs(d) > 1e-12 else 0.0, -0.5, 0.5))
    return (i + off) * df, spec[i] / floor, view


def autocorr_pitch(x, fs, flo, fhi):
    n = len(x)
    xw = (x - np.mean(x)) * np.hanning(n)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    spec = np.fft.rfft(xw, nfft)
    ac = np.fft.irfft(spec * np.conj(spec))[:n]
    if ac[0] <= 0:
        return None, 0.0
    ac /= ac[0]
    lag_min = max(2, int(fs / fhi)); lag_max = min(n - 2, int(fs / flo))
    if lag_max <= lag_min:
        return None, 0.0
    seg = ac[lag_min:lag_max]
    best = float(np.max(seg))
    peaks = np.where((seg[1:-1] >= seg[:-2]) & (seg[1:-1] >= seg[2:]) &
                     (seg[1:-1] >= 0.85 * best))[0] + 1
    i = int(peaks[0]) if len(peaks) else int(np.argmax(seg))
    k = i + lag_min
    y0, y1, y2 = ac[k - 1], ac[k], ac[k + 1]
    d = y0 - 2 * y1 + y2
    k += float(np.clip(0.5 * (y0 - y2) / d if abs(d) > 1e-12 else 0.0, -0.5, 0.5))
    return fs / k, float(seg[i])


def tone_in_window(x, fs, fmin, fmax, min_prom=3.0):
    """One sliding window -> (f or None, prominence, spectrum view)."""
    x = x - np.mean(x)
    f, prom, view = spectral_peak(x, fs, fmin, fmax)
    if f is None or prom < min_prom:
        return None, prom, view
    xb = bandpass_fft(x, fs, f * 0.9, f * 1.1)
    f_ac, clarity = autocorr_pitch(xb, fs, f * 0.85, f * 1.15)
    if f_ac is None or clarity < 0.5 or abs(f_ac - f) / f > 0.03:
        return None, prom, view
    return f, prom, view


class ToneDetector:
    """One reading per pluck. On an abrupt onset it starts TRACKING: every window's
    tone estimate is collected until the note dies away (or MAX_TRACK_S elapses).
    The reading is the frequency cluster that persisted for the most windows, with
    ties going to the later cluster, so short-lived attack components lose to the
    sustained segment tone. While tracking, nothing else can fire; a genuinely new
    onset (a sharp step well above the current level) ends the track early and
    starts the next one, so quick successive plucks still each count once."""

    ONSET_RATIO = 4.0
    ONSET_SHARPNESS = 6.0
    SETTLE_HOPS = 2
    MAX_TRACK_S = 0.7
    END_QUIET_HOPS = 3            # windows without a tone -> note has ended
    MIN_TONE_FRACTION = 0.15
    MIN_CLUSTER = 2

    def __init__(self, fs, fmin, fmax):
        self.fs, self.fmin, self.fmax = fs, fmin, fmax
        self.energy = []
        self.past = []
        self.tracking = False
        self.age = 10 ** 6
        self.track = []            # (age, f, amp, frac) while tracking
        self.no_tone = 0
        self.sharp = 0.0
        self.pre_x = None
        self.last_view = None
        self.accept_view = None
        self.last_prom = 0.0
        self.last_reason = "idle"
        self.hop_s = 0.04
        self.expected_f = None     # prior from the wheel (side mean) for the spoke being read

    # ---------------------------------------------------------------- helpers
    def _onset(self, xb, e):
        prev = self.energy[-8:-1]
        if len(prev) < 3:
            return False
        if self.tracking:
            # a new pluck on top of a ringing one must stand well clear of the current level
            return e > self.ONSET_RATIO * self.energy[-1] and self.age > self.SETTLE_HOPS + 1
        rising = e > 1.5 * self.energy[-1]
        return rising and e > self.ONSET_RATIO * float(np.median(prev))

    def _sharpness(self, xb):
        sub = int(0.01 * self.fs)
        E = np.array([np.mean(xb[i:i + sub] ** 2) for i in range(0, len(xb) - sub, sub)]) + 1e-12
        floor = float(np.median(E))
        return max((E[i] / max(E[i - 1], floor) for i in range(1, len(E)) if E[i] > 4 * floor),
                   default=0.0)

    def _finish(self):
        """Pick the reading from the collected track. Returns Hz or None."""
        self.tracking = False
        pts = [p for p in self.track if p[0] >= self.SETTLE_HOPS]
        if len(pts) < self.MIN_CLUSTER:               # very short note: use the earlier windows too
            pts = [p for p in self.track if p[0] >= 1]
        self.track = []
        if self.sharp < self.ONSET_SHARPNESS:
            self.last_reason = f"onset too soft ({self.sharp:.1f}x)"
            return None
        if len(pts) < self.MIN_CLUSTER:
            self.last_reason = "note too short"
            return None
        # cluster by frequency (3 %)
        clusters = []
        for age, f, amp, frac in sorted(pts, key=lambda p: p[1]):
            if clusters and abs(f - clusters[-1][-1][1]) / f < 0.03:
                clusters[-1].append((age, f, amp, frac))
            else:
                clusters.append([(age, f, amp, frac)])
        def cluster_score(c):
            fc = float(np.median([p[1] for p in c]))
            score = max(p[2] for p in c) * np.sqrt(len(c))  # strength x (some) persistence
            if self.expected_f:
                dev = abs(fc - self.expected_f) / self.expected_f
                if dev > 0.20:                               # far from the wheel's prior: needs ~3x the evidence
                    score *= 0.3
            return (score, np.mean([p[0] for p in c]))
        clusters.sort(key=cluster_score, reverse=True)
        best = clusters[0]
        if len(best) < self.MIN_CLUSTER:
            self.last_reason = "no persistent tone"
            return None
        f = float(np.median([p[1] for p in best]))
        if max(p[3] for p in best) < self.MIN_TONE_FRACTION:
            self.last_reason = "tone weak vs background"
            return None
        amps = [p[2] for p in sorted(best)]
        if len(amps) >= 3 and amps[-1] > 1.3 * amps[0]:
            self.last_reason = "tone growing, not a pluck"
            return None
        if self.pre_x is not None:
            pt = bandpass_fft(self.pre_x, self.fs, f * 0.9, f * 1.1)
            if max(amps) < 3.0 * float(np.sqrt(np.mean(pt * pt))):
                self.last_reason = "tone present before onset (hum)"
                return None
        self.last_reason = f"pluck ({len(best)}/{len(pts)} windows)"
        return f

    # ------------------------------------------------------------------ main
    def feed(self, x):
        x = x - np.mean(x)
        xb = bandpass_fft(x, self.fs, self.fmin * 0.9, self.fmax * 1.1)
        e = float(np.mean(xb * xb)) + 1e-12
        result = None

        if self._onset(xb, e):
            if self.tracking:
                result = self._finish()          # previous note ends now
            self.tracking = True
            self.age = 0
            self.track = []
            self.no_tone = 0
            self.sharp = self._sharpness(xb)
            self.pre_x = self.past[-5] if len(self.past) >= 5 else (self.past[0] if self.past else None)
        else:
            self.age += 1
        self.energy = (self.energy + [e])[-12:]
        self.past = (self.past + [x])[-6:]

        f, prom, view = tone_in_window(x, self.fs, self.fmin, self.fmax)
        self.last_view, self.last_prom = view, prom

        if not self.tracking:
            self.last_reason = "idle" if f is None else "tone (no onset)"
            return result

        if f is None:
            self.no_tone += 1
        else:
            self.no_tone = 0
            xt = bandpass_fft(x, self.fs, f * 0.9, f * 1.1)
            xa = bandpass_fft(x, self.fs, 100.0, self.fs / 2 - 100.0)
            tail = len(x) // 3
            frac = float(np.mean(xt[tail:] ** 2) / (np.mean(xa[tail:] ** 2) + 1e-12))
            if frac < 0.05:                      # tone buried in noise: the note has effectively ended
                self.no_tone += 1
            else:
                self.track.append((self.age, f, float(np.sqrt(np.mean(xt * xt))), frac))
                if self.age >= self.SETTLE_HOPS:
                    self.accept_view = view
        self.last_reason = f"tracking ({self.age * self.hop_s:.2f}s)"

        ended = (self.no_tone >= self.END_QUIET_HOPS and self.age > self.SETTLE_HOPS + self.END_QUIET_HOPS)
        if ended or self.age * self.hop_s >= self.MAX_TRACK_S:
            result = self._finish() if result is None else result
        return result


# ---------------------------------------------------------------- session
class Wheel:
    def __init__(self, n, first_side, ratio, path):
        self.n = n
        self.sides = [first_side if i % 2 == 0 else ("R" if first_side == "L" else "L")
                      for i in range(n)]
        self.ratio = ratio          # expected T_R / T_L from flange geometry
        self.t = [None] * n         # kgf per spoke, 1-based in UI
        self.f = [None] * n         # Hz per spoke
        self.spec = [None] * n      # stored spectrum at 1 Hz bins, 0-2000 Hz, normalised
        self.cur = 0
        self.gain = 40.0            # kgf per full turn, refined by calibration
        self.gain_samples = []
        self.calib = None           # pending calibration (spoke idx, kgf before)
        self.path = path
        self.undo = []
        self.length_m = None
        self.pass_log = []          # [(mean_L, mean_R, turns_L, turns_R)] per completed pass
        self.suspect = set()        # spoke indices whose last reading was far from the side mean
        if path and os.path.exists(path):
            self.load()

    def side_name(self, i):
        return "disc" if self.sides[i] == "L" else "non-disc"

    def record(self, kgf, f=None, view=None):
        self.undo.append((self.cur, self.t[self.cur], self.f[self.cur], self.spec[self.cur]))
        self.t[self.cur] = kgf
        self.f[self.cur] = f
        if view is not None:
            vf, vm = view
            grid = np.arange(0, 2000)
            sp = np.interp(grid, vf, vm)
            self.spec[self.cur] = [round(float(v), 4) for v in sp / (np.max(sp) + 1e-12)]
        self.save()

    def undo_last(self):
        if self.undo:
            i, v, f, sp = self.undo.pop()
            self.t[i], self.f[i], self.spec[i] = v, f, sp
            self.cur = i
            self.save()

    def side_stats(self, side):
        vals = [v for v, s in zip(self.t, self.sides) if s == side and v is not None]
        return (float(np.mean(vals)), float(np.std(vals)), len(vals)) if vals else (None, None, 0)

    def recommendations(self):
        """Turns to bring each measured spoke to its side's mean. +ve = tighten."""
        out = []
        for i, (v, s) in enumerate(zip(self.t, self.sides)):
            if v is None:
                continue
            m, _, _ = self.side_stats(s)
            d = m - v
            out.append((i, v, d, d / self.gain))
        out.sort(key=lambda r: -abs(r[2]))
        return out

    def save(self):
        """Safe save: write new file -> read it back and verify -> rotate old to .bak -> rename new in.
        The live session file is never opened for writing, so an interruption can't truncate it."""
        if not self.path:
            return False
        data = {"n": self.n, "sides": self.sides, "t": self.t, "gain": self.gain,
                "gain_samples": self.gain_samples, "ratio": self.ratio,
                "length_m": self.length_m, "f": self.f, "spec": self.spec,
                "pass_log": self.pass_log}
        tmp = self.path + ".tmp"
        try:
            with open(tmp, "w") as fh:
                json.dump(data, fh, indent=1)
                fh.flush(); os.fsync(fh.fileno())
            # verify by reading it back
            with open(tmp) as fh:
                check = json.load(fh)
            ok = (check.get("n") == self.n and len(check.get("t", [])) == self.n
                  and len(check.get("sides", [])) == self.n and check.get("t") == self.t)
            if not ok:
                raise ValueError("verification mismatch")
        except Exception as ex:
            print(f"save aborted, previous session kept: {ex}")
            try:
                os.remove(tmp)
            except OSError:
                pass
            return False
        if os.path.exists(self.path) and os.path.getsize(self.path) > 0:
            os.replace(self.path, self.path + ".bak")
        os.replace(tmp, self.path)
        return True

    def load(self):
        """Load the session; fall back to the .bak if the main file is empty or corrupt."""
        for candidate in (self.path, self.path + ".bak"):
            if not os.path.exists(candidate):
                continue
            try:
                d = json.load(open(candidate))
            except (json.JSONDecodeError, OSError):
                print(f"warning: {candidate} is empty or corrupt, skipping it")
                continue
            if d.get("n") != self.n:
                continue
            self.sides, self.t = d["sides"], d["t"]
            self.gain = d.get("gain", self.gain)
            self.gain_samples = d.get("gain_samples", [])
            self.length_m = d.get("length_m")
            self.f = d.get("f", [None] * self.n)
            self.spec = d.get("spec", [None] * self.n)
            self.pass_log = d.get("pass_log", [])
            if candidate != self.path:
                print(f"recovered session from {candidate}")
            return
        print("starting a fresh session")


def fmt_turns(turns):
    """Express turns as the nearest eighth, e.g. '+1/4', '-3/8', '+1 1/8'."""
    e = int(round(turns * 8))
    if e == 0:
        return "0"
    sign = "+" if e > 0 else "-"
    e = abs(e)
    whole, frac = divmod(e, 8)
    fr = {0: "", 1: "1/8", 2: "1/4", 3: "3/8", 4: "1/2", 5: "5/8", 6: "3/4", 7: "7/8"}[frac]
    return sign + (f"{whole} " if whole else "") + fr if fr else f"{sign}{whole}"


# ------------------------------------------------------------- user config
CONFIG_PATH = os.path.expanduser("~/.spoke_tension_config.json")


def load_config():
    try:
        return json.load(open(CONFIG_PATH))
    except Exception:
        return {}


def save_config(**kw):
    cfg = load_config(); cfg.update(kw)
    try:
        tmp = CONFIG_PATH + ".tmp"
        json.dump(cfg, open(tmp, "w"), indent=1); os.replace(tmp, CONFIG_PATH)
    except Exception as ex:
        print("could not write config:", ex)


# ------------------------------------------------------------- wav helpers
def write_wav(path, x, fs):
    import wave
    pcm = (np.clip(x, -1, 1) * 32767).astype("<i2")
    with wave.open(path, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(fs); w.writeframes(pcm.tobytes())


def read_wav(path):
    import wave
    with wave.open(path, "rb") as w:
        fs, n, ch, sw = w.getframerate(), w.getnframes(), w.getnchannels(), w.getsampwidth()
        raw = w.readframes(n)
    x = np.frombuffer(raw, dtype={1: np.int8, 2: "<i2", 4: "<i4"}[sw]).astype(np.float64)
    if ch > 1:
        x = x.reshape(-1, ch)[:, 0]
    return x / float(2 ** (8 * sw - 1)), fs


def analyse_file(path, length_m, mu, fmin, fmax):
    x, fs = read_wav(path)
    det = ToneDetector(fs, fmin, fmax)
    hop = int(0.04 * fs); win = int(TONE_WIN_S * fs)
    hits = []
    for a in range(0, len(x) - win, hop):
        f = det.feed(x[a:a + win])
        if f is not None:
            hits.append((a / fs, f))
    print(f"{path}: fs={fs}, {len(hits)} pluck(s) detected")
    for t, f in hits:
        print(f"  t={t:.2f}s  {f:.1f} Hz -> {freq_to_kgf(f, length_m, mu):.1f} kgf")


# ---------------------------------------------------------------- devices
def choose_device(preselected):
    devices = sd.query_devices()
    inputs = [(i, d) for i, d in enumerate(devices) if d["max_input_channels"] > 0]
    if not inputs:
        sys.exit("No input devices found.")
    if preselected is not None:
        if preselected.isdigit():
            return int(preselected)
        def plain(txt):
            import unicodedata
            txt = txt.replace("Ø", "O").replace("ø", "o")
            txt = unicodedata.normalize("NFKD", txt)
            return "".join(c for c in txt if not unicodedata.combining(c)).lower()
        matches = [i for i, d in inputs if plain(preselected) in plain(d["name"])]
        if not matches:
            sys.exit(f"no input device matching '{preselected}'. Available: "
                     + ", ".join(d["name"] for _, d in inputs))
        if len(matches) > 1:
            print(f"note: several devices match '{preselected}', using the first: "
                  f"{devices[matches[0]]['name']}")
        return matches[0]
    default_in = sd.default.device[0]
    print("\nInput devices:")
    for i, d in inputs:
        print(f"  [{i}] {d['name']}  ({int(d['default_samplerate'])} Hz)"
              + (" (default)" if i == default_in else ""))
    while True:
        s = input(f"\nSelect device [{default_in}]: ").strip()
        if s == "":
            return default_in
        if s.isdigit() and int(s) in [i for i, _ in inputs]:
            return int(s)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="Spoke tension by pluck frequency, with wheel map.")
    ap.add_argument("--length", type=float, default=186.0, help="vibrating length mm (default 186)")
    ap.add_argument("--full-length", type=float, default=260.0,
                    help="free length nipple-to-flange when the crossing partner is held off (mm)")
    ap.add_argument("--fret-length", type=float, default=186.0,
                    help="physical crossing-to-nipple distance (mm), used by the fret calibration method")
    ap.add_argument("--diameter", type=float, default=1.8, help="spoke diameter mm (default 1.8)")
    ap.add_argument("--tmin", type=float, default=20.0)
    ap.add_argument("--tmax", type=float, default=140.0)
    ap.add_argument("--spokes", type=int, default=28)
    ap.add_argument("--first-side", choices=["L", "R"], default="L",
                    help="side of spoke 1 (the hole next to the valve); L = disc side")
    ap.add_argument("--flange-left", type=float, default=22.2, help="flange distance left mm")
    ap.add_argument("--flange-right", type=float, default=31.3, help="flange distance right mm")
    ap.add_argument("--session", type=str, default="wheel_session.json")
    ap.add_argument("--target", type=float, default=None,
                    help="final disc-side tension in kgf (e.g. 115); the non-disc target follows from the ratio")
    ap.add_argument("--max-pass", type=float, default=0.25,
                    help="largest uniform adjustment per pass, in turns (default 1/4)")
    ap.add_argument("--device", type=str, default=None,
                    help="input device index, or part of its name, e.g. --device rode")
    ap.add_argument("--rate", type=int, default=None)
    ap.add_argument("--gain", type=float, default=1.0, help="software mic gain")
    ap.add_argument("--record", type=str, default=None, help="save each accepted pluck as wav")
    ap.add_argument("--file", type=str, default=None, help="analyse a wav offline")
    args = ap.parse_args()

    length_m = args.length / 1000.0
    mu = mass_per_metre(args.diameter)
    fmin, fmax = kgf_to_freq(args.tmin, length_m, mu), kgf_to_freq(args.tmax, length_m, mu)

    if args.file:
        analyse_file(args.file, length_m, mu, fmin, fmax)
        return
    geo = {"L": length_m, "Lfull": args.full_length / 1000.0, "lcal": None, "lcal_f": []}
    SEGMENT_RATIO[0] = geo["Lfull"] / geo["L"]

    dev = choose_device(args.device)
    info = sd.query_devices(dev)
    fs = args.rate or int(info["default_samplerate"])
    if args.record:
        os.makedirs(args.record, exist_ok=True)

    ratio = args.flange_left / args.flange_right   # T_R / T_L (small-angle approximation)
    wheel = Wheel(args.spokes, args.first_side, ratio, args.session)
    cfg = load_config()
    length_src = "default"
    if args.length != 186.0:
        length_src = "--length"
    elif wheel.length_m:
        geo["L"] = wheel.length_m; length_src = "session"
    elif cfg.get("length_m"):
        geo["L"] = cfg["length_m"]; wheel.length_m = geo["L"]; length_src = "config"
    lo_ok, hi_ok = args.fret_length / 1000.0 * 0.95, geo["Lfull"] * 1.05
    if not (lo_ok <= geo["L"] <= hi_ok):
        print(f"\n*** stored length {geo['L']*1000:.0f} mm is impossible (must be {lo_ok*1000:.0f}-{hi_ok*1000:.0f} mm); "
              f"ignoring it. Pass --length or calibrate. ***")
        geo["L"] = args.length / 1000.0; wheel.length_m = None; length_src = "default"
    SEGMENT_RATIO[0] = geo["Lfull"] / geo["L"]
    if length_src == "--length":
        wheel.length_m = geo["L"]; save_config(length_m=geo["L"], full_length_mm=args.full_length)
    for i_, fi in enumerate(wheel.f):          # keep kgf consistent with the length in use
        if fi:
            wheel.t[i_] = freq_to_kgf(fi, geo["L"], mu)
    # detector band must cover the full-length tone too (lower by Lfull/L)
    det = ToneDetector(fs, fmin, fmax)
    band_normal = (fmin, fmax)
    band_wide = (fmin * geo["L"] / geo["Lfull"] * 0.9, fmax)

    print(f"\nDevice {dev}: {info['name']} @ {fs} Hz")
    if length_src == "default":
        print("\n" + "*" * 70)
        print("*** VIBRATING LENGTH IS THE 186 mm DEFAULT: readings will be LOW.          ***")
        print("*** Press l in the window to calibrate (two plucks), or pass --length NNN  ***")
        print("*" * 70)
    else:
        print(f"Vibrating length {geo['L']*1000:.0f} mm ({length_src})")
    print(f"d = {args.diameter} mm; band {fmin:.0f}-{fmax:.0f} Hz ({args.tmin:.0f}-{args.tmax:.0f} kgf)")
    print(f"Expected non-disc/disc tension ratio from flange geometry: {ratio:.2f}")
    print(f"Session: {args.session}   (starts PAUSED: press space in the window to listen)\n")

    buf_len = int(fs * 2.0)
    buf = np.zeros(buf_len, dtype=np.float32)
    st = {"armed": False, "msg": "paused - press space in the window to start", "digits": "",
          "nrec": 0, "last": None, "last_t": 0.0, "plan": (0.0, 0.0)}

    def callback(indata, frames, t, status):
        nonlocal buf
        buf = np.roll(buf, -frames)
        buf[-frames:] = indata[:, 0] * args.gain

    stream = sd.InputStream(device=dev, channels=1, samplerate=fs, blocksize=1024, callback=callback)
    # the stream only runs while listening; paused = no capture, no analysis

    def set_listening(on):
        nonlocal buf
        if on and not stream.active:
            buf[:] = 0.0                      # start from silence, not stale audio
            det.__init__(fs, det.fmin, det.fmax)
            stream.start()
        elif not on and stream.active:
            stream.abort()
        st["armed"] = on

    shutting_down = {"done": False}

    def shutdown(*_):
        """Idempotent, fast exit: abort audio (no drain), save once, close the window."""
        if shutting_down["done"]:
            return
        shutting_down["done"] = True
        try:
            if stream.active:
                stream.abort()       # immediate, unlike stop() which waits for buffers
            stream.close()
        except Exception:
            pass
        try:
            wheel.save()
            print(f"\nsession saved to {wheel.path}")
        except Exception as ex:
            print("save failed:", ex)
        try:
            plt.close("all")
        except Exception:
            pass

    # Ctrl-C in the terminal -> same clean path as pressing q
    signal.signal(signal.SIGINT, shutdown)
    signal.signal(signal.SIGTERM, shutdown)

    # ---------------------------------------------------------------- figure
    plt.rcParams["toolbar"] = "None"
    fig = plt.figure(figsize=(15, 9))
    fig.canvas.manager.set_window_title("Spoke tension")
    outer = fig.add_gridspec(1, 2, width_ratios=[3.2, 1.5], wspace=0.08,
                             left=0.045, right=0.985, top=0.95, bottom=0.07)
    left = outer[0].subgridspec(2, 2, height_ratios=[2.6, 1.3], width_ratios=[1, 2.2], hspace=0.3, wspace=0.15)
    right = outer[1].subgridspec(3, 1, height_ratios=[1.0, 3.2, 1.3], hspace=0.12)

    # ---- wheel map (big)
    ax_map = fig.add_subplot(left[0, :])
    ax_map.set_title("Wheel map   blue = disc side   orange = non-disc side   gold = current spoke   hatched = suspect, re-pluck", fontsize=11)
    idx = np.arange(1, wheel.n + 1)
    colors = ["tab:blue" if sd_ == "L" else "tab:orange" for sd_ in wheel.sides]
    bars = ax_map.bar(idx, [0] * wheel.n, color=colors, alpha=0.85)
    ax_map.set_xticks(idx); ax_map.set_xlim(0.3, wheel.n + 0.7)
    ax_map.set_ylim(0, args.tmax); ax_map.set_ylabel("kgf"); ax_map.set_xlabel("spoke (1 = beside valve, clockwise)")
    ax_map.grid(axis="y", alpha=0.3)
    mean_L = ax_map.axhline(0, color="tab:blue", lw=1.2, ls="--", visible=False)
    mean_R = ax_map.axhline(0, color="tab:orange", lw=1.2, ls="--", visible=False)
    cur_marker = ax_map.axvspan(0.5, 1.5, color="gold", alpha=0.35)
    if args.target:
        ax_map.axhline(args.target, color="tab:blue", lw=1.5, alpha=0.6)
        ax_map.axhline(args.target * ratio, color="tab:orange", lw=1.5, alpha=0.6)
        ax_map.text(wheel.n + 0.6, args.target, f" target {args.target:.0f}", color="tab:blue", fontsize=8, va="center")
        ax_map.text(wheel.n + 0.6, args.target * ratio, f" target {args.target*ratio:.0f}", color="tab:orange", fontsize=8, va="center")
    bar_labels = [ax_map.text(i, 0, "", ha="center", va="bottom", fontsize=8) for i in idx]

    # ---- small live plots
    ax_wave = fig.add_subplot(left[1, 0]); ax_wave.set_title("mic (off while paused)", fontsize=9)
    nwin = int(fs * 0.5)
    (wave_line,) = ax_wave.plot(np.linspace(0, 500, nwin), np.zeros(nwin), lw=0.5)
    ax_wave.set_xlim(0, 500); ax_wave.set_yticks([]); ax_wave.tick_params(labelsize=7)
    ax_spec = fig.add_subplot(left[1, 1])
    ax_spec.set_title("spectrum   blue = live   red = stored reading of selected spoke", fontsize=9)
    ax_spec.set_xlim(0, 1500); ax_spec.set_yticks([]); ax_spec.tick_params(labelsize=7)
    ax_spec.axvspan(fmin, fmax, color="tab:green", alpha=0.15)
    (spec_line,) = ax_spec.plot([], [], lw=0.7, color="tab:blue", alpha=0.7)
    (spec_stored,) = ax_spec.plot([], [], lw=1.0, color="tab:red")
    spec_fline = ax_spec.axvline(0, color="tab:red", lw=1, ls="--", visible=False)
    spec_label = ax_spec.text(0.99, 0.95, "", transform=ax_spec.transAxes, ha="right", va="top",
                              fontsize=8, color="tab:red")

    def panel(cell, title):
        ax = fig.add_subplot(cell)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("0.75")
        ax.set_title(title, fontsize=9, loc="left", color="0.35")
        return ax

    ax_read = panel(right[0], "last reading")
    read_big = ax_read.text(0.5, 0.56, "", ha="center", va="center", fontsize=30, weight="bold", family="monospace")
    read_sub = ax_read.text(0.5, 0.15, "", ha="center", va="center", fontsize=10, color="0.4", family="monospace")
    state_txt = ax_read.text(0.02, 0.92, "", ha="left", va="top", fontsize=10, weight="bold")

    ax_info = panel(right[1], "wheel")
    info_t = ax_info.text(0.03, 0.97, "", va="top", fontsize=9.2, family="monospace")

    ax_keys = panel(right[2], "keys")
    ax_keys.text(0.03, 0.95, "\n".join([
        "space      start / stop listening",
        "enter  ->  next spoke     <-  previous",
        "N enter    jump to spoke N",
        "c          calibrate turn gain (+1/4 turn)",
        "l          calibrate length (guided: free or fret method)",
        "p          log a uniform pass (then re-read all)",
        "u  undo    s  save    q / esc / Ctrl-C  quit (saves)"]),
        va="top", fontsize=9, family="monospace")

    status = fig.text(0.045, 0.012, "", fontsize=8.5, color="0.35")

    # calibration overlay (hidden unless a calibration is running)
    overlay = ax_map.text(0.5, 0.55, "", transform=ax_map.transAxes, ha="center", va="center",
                          fontsize=13, family="monospace", visible=False, zorder=20,
                          bbox=dict(boxstyle="round,pad=1.0", fc="#fff8dc", ec="tab:red", lw=2))

    def set_overlay(lines):
        if lines:
            overlay.set_text("\n".join(lines)); overlay.set_visible(True)
        else:
            overlay.set_visible(False)


    def propose_length(L_eff, how_lines):
        """Bounds-check a calibrated length and ask for confirmation before applying it."""
        lo = args.fret_length / 1000.0 * 0.95
        hi = geo["Lfull"] * 1.05
        if not (lo <= L_eff <= hi):
            geo["lcal_msg"] = how_lines + ["",
                f"RESULT {L_eff*1000:.0f} mm IS IMPOSSIBLE: it must lie between the",
                f"crossing-to-nipple distance ({lo*1000:.0f} mm) and the whole",
                f"spoke ({hi*1000:.0f} mm). One of the notes was not the spoke.",
                "", "No change made.  Press l to try again."]
            geo["lcal"] = "done"
            return
        ex = next(((i, fi) for i, fi in enumerate(wheel.f) if fi), None)
        example = []
        if ex:
            i, fi = ex
            example = [f"e.g. spoke {i+1} ({fi:.0f} Hz): {freq_to_kgf(fi, geo['L'], mu):.0f} kgf -> "
                       f"{freq_to_kgf(fi, L_eff, mu):.0f} kgf"]
        geo["lcal_pending"] = L_eff
        geo["lcal_msg"] = how_lines + ["",
            f"Proposed working length: {L_eff*1000:.0f} mm  (currently {geo['L']*1000:.0f} mm)"] + example + [
            "", "Apply this?   y = yes, save and re-convert all readings",
            "              n = no, keep the current length"]
        geo["lcal"] = "confirm"

    def apply_length(L_eff):
        old_L = geo["L"]
        geo["L"] = L_eff; wheel.length_m = L_eff
        SEGMENT_RATIO[0] = geo["Lfull"] / geo["L"]
        for i_, fi in enumerate(wheel.f):
            if fi:
                wheel.t[i_] = freq_to_kgf(fi, L_eff, mu)
        wheel.save()
        save_config(length_m=L_eff, full_length_mm=args.full_length)
        print(f"length set to {L_eff*1000:.0f} mm (was {old_L*1000:.0f} mm); readings re-converted")

    def cal_instructions():
        step = geo["lcal"]
        head = "LENGTH CALIBRATION          (esc to cancel)"
        if step == "choose":
            return [head, "",
                    "Choose a method:",
                    "",
                    "  1   FREE method - hold the crossing partner clear so the",
                    f"      whole spoke rings ({geo['Lfull']*1000:.0f} mm). Easy at low tension.",
                    "",
                    "  2   FRET method - pinch the crossing and brace your hand on",
                    "      the rim so the crossing cannot move, then pluck",
                    f"      ({args.fret_length:.0f} mm crossing-to-nipple). Fiddly at tension.",
                    "",
                    "  3   KNOWN TENSION - you have a tensiometer reading for the",
                    "      current spoke. Type it in; one pluck sets the length.",
                    "      Most reliable.",
                    "",
                    "Press 1, 2 or 3"]
        if step == "known_enter":
            return [head, "",
                    "KNOWN TENSION  -  STEP 1 of 2",
                    "",
                    f"Type the measured tension of spoke {wheel.cur+1} in kgf,",
                    "then press enter.",
                    "",
                    f"    > {geo['lcal_typed']}_"]
        if step == "known_pluck":
            return [head, "",
                    f"Known tension: {geo['lcal_T']:.0f} kgf",
                    "",
                    "KNOWN TENSION  -  STEP 2 of 2",
                    "",
                    f"Pluck spoke {wheel.cur+1} normally.",
                    "", "Listening..."]
        if step == "fret_normal":
            return [head, "",
                    "FRET METHOD  -  STEP 1 of 2:  normal pluck",
                    "",
                    "Pluck the spoke as you normally would,",
                    "touching nothing.",
                    "", "Listening for the normal note..."]
        if step == "fret_clamped":
            return [head, "",
                    f"Normal note: {geo['lcal_f'][0]:.0f} Hz   OK",
                    "",
                    "FRET METHOD  -  STEP 2 of 2:  clamped pluck",
                    "",
                    "Press the spoke HARD against its crossing partner,",
                    "exactly at the crossing, with a fingernail or hex key.",
                    "Hold both spokes still there and pluck the SAME spoke",
                    "between the crossing and the rim.",
                    "", "Listening for the clamped note..."]
        if step == "free":
            return [head, "",
                    "STEP 1 of 2:  measure the FULL spoke",
                    "",
                    "Put a finger on the crossing partner and pull it",
                    "clear of this spoke so they do not touch.",
                    f"Pluck the spoke (full free length {geo['Lfull']*1000:.0f} mm).",
                    "", "Listening for the full-length note..."]
        if step == "crossed":
            return [head, "",
                    f"Full-length note: {geo['lcal_f'][0]:.0f} Hz   OK",
                    "",
                    "STEP 2 of 2:  measure the CROSSED spoke",
                    "",
                    "Let go.  Pluck the SAME spoke normally,",
                    "with the crossing touching as laced.",
                    "", "Listening for the crossed note..."]
        if step == "confirm":
            return [head, ""] + geo["lcal_msg"]
        if step == "done":
            return [head, ""] + geo["lcal_msg"] + ["", "Press any key to continue"]
        return []

    # ------------------------------------------------------------- key input
    def on_key(ev):
        k = ev.key
        if geo["lcal"] == "confirm":
            if k == "y":
                apply_length(geo["lcal_pending"])
                geo["lcal_msg"] = [f"Working length set to {geo['L']*1000:.0f} mm.",
                                   "Saved to session and config; all readings re-converted."]
                geo["lcal"] = "done"; set_overlay(cal_instructions())
            elif k in ("n", "escape"):
                geo["lcal"] = None; set_overlay(None); set_listening(geo.get("armed_before", False))
                st["msg"] = "calibration discarded, length unchanged"
            elif k == "q":
                shutdown()
            return
        if geo["lcal"] == "done":                     # dismiss result screen
            geo["lcal"] = None; set_overlay(None)
            set_listening(geo.get("armed_before", False))
            return
        if geo["lcal"] == "known_enter":
            if k and (k.isdigit() or k == "."):
                geo["lcal_typed"] += k
            elif k == "backspace":
                geo["lcal_typed"] = geo["lcal_typed"][:-1]
            elif k == "enter":
                try:
                    T = float(geo["lcal_typed"])
                    if not (10 <= T <= 200):
                        raise ValueError
                    geo["lcal_T"] = T; geo["lcal"] = "known_pluck"
                except ValueError:
                    geo["lcal_typed"] = ""
            elif k in ("escape", "l"):
                geo["lcal"] = None; set_overlay(None); set_listening(geo.get("armed_before", False))
                return
            elif k == "q":
                shutdown(); return
            set_overlay(cal_instructions())
            return
        if geo["lcal"] == "choose":
            if k == "3":
                geo["lcal"] = "known_enter"; geo["lcal_typed"] = ""; set_overlay(cal_instructions())
            elif k == "1":
                geo["lcal"] = "free"; det.fmin, det.fmax = band_wide; set_overlay(cal_instructions())
            elif k == "2":
                geo["lcal"] = "fret_normal"; set_overlay(cal_instructions())
            elif k in ("escape", "l"):
                geo["lcal"] = None; set_overlay(None); set_listening(geo.get("armed_before", False))
            elif k == "q":
                shutdown()
            return
        if geo["lcal"] is not None:                   # calibrating: only esc / q do anything
            if k in ("escape", "l"):
                geo["lcal"] = None; det.fmin, det.fmax = band_normal; set_overlay(None)
                set_listening(geo.get("armed_before", False)); st["msg"] = "length calibration cancelled"
            elif k == "q":
                shutdown()
            return
        if k == " ":
            set_listening(not st["armed"])
            st["msg"] = "listening" if st["armed"] else "paused"
        elif k in ("enter", "right"):
            if st["digits"]:
                n = int(st["digits"]); st["digits"] = ""
                if 1 <= n <= wheel.n:
                    wheel.cur = n - 1
            else:
                wheel.cur = (wheel.cur + 1) % wheel.n
        elif k == "left":
            wheel.cur = (wheel.cur - 1) % wheel.n
        elif k and k.isdigit():
            st["digits"] += k
        elif k == "u":
            wheel.undo_last(); st["msg"] = "undid last reading"
        elif k == "s":
            wheel.save(); st["msg"] = f"saved {wheel.path}"
        elif k == "c":
            v = wheel.t[wheel.cur]
            if wheel.calib is None:
                if v is None:
                    st["msg"] = "measure this spoke first, then press c"
                else:
                    wheel.calib = (wheel.cur, v)
                    st["msg"] = f"GAIN CAL: turn spoke {wheel.cur+1} +1/4 turn, re-pluck, press c"
            else:
                i, before = wheel.calib
                if i == wheel.cur and v is not None and v != before:
                    g = (v - before) / 0.25
                    if g > 0:
                        wheel.gain_samples.append(g)
                        wheel.gain = float(np.mean(wheel.gain_samples))
                        st["msg"] = f"gain now {wheel.gain:.0f} kgf/turn ({len(wheel.gain_samples)} samples)"
                    else:
                        st["msg"] = "tension went down? calibration discarded"
                else:
                    st["msg"] = "calibration cancelled"
                wheel.calib = None
                wheel.save()
        elif k == "l":
            geo["lcal"] = "choose"; geo["lcal_f"] = []
            geo["armed_before"] = st["armed"]
            set_listening(True)
            st["msg"] = "length calibration running - spoke readings are paused"
            set_overlay(cal_instructions())
        elif k == "p":
            mL, _, nL = wheel.side_stats("L"); mR, _, nR = wheel.side_stats("R")
            if nL + nR < wheel.n:
                st["msg"] = "read every spoke before logging a pass"
            else:
                tl, tr = st["plan"]
                wheel.pass_log.append([mL, mR, tl, tr]); wheel.save()
                st["msg"] = (f"pass logged: disc {fmt_turns(tl)}, non-disc {fmt_turns(tr)} turn applied. "
                             f"Now re-read all 28 spokes")
        elif k in ("q", "escape"):
            shutdown()

    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("close_event", shutdown)
    # a periodic no-op timer lets Python service SIGINT while the GUI loop is blocking
    wake = fig.canvas.new_timer(interval=200)
    wake.add_callback(lambda: None)
    wake.start()

    # ------------------------------------------------------------- animation
    def update(_):
        now = time.time()
        f = None
        m_cur, _, n_cur = wheel.side_stats(wheel.sides[wheel.cur])
        if st["armed"]:
            x = buf[-nwin:].astype(np.float64)
            wave_line.set_ydata(x)
            amp = max(np.max(np.abs(x)) * 1.2, 1e-3); ax_wave.set_ylim(-amp, amp)
            win = buf[-int(TONE_WIN_S * fs):].astype(np.float64)
            det.expected_f = kgf_to_freq(m_cur, geo["L"], mu) if (m_cur and n_cur >= 3 and geo["lcal"] is None) else None
            f = det.feed(win)
            if det.last_view is not None:
                vf, vm = det.last_view
                spec_line.set_data(vf, vm / (np.max(vm) + 1e-12)); ax_spec.set_ylim(0, 1.05)
        else:
            wave_line.set_ydata(np.zeros(nwin)); ax_wave.set_ylim(-1e-3, 1e-3)
            spec_line.set_data([], [])

        if geo["lcal"] is not None:
            # ---- calibration owns every pluck; nothing is routed to spokes
            if geo["lcal"] == "choose":
                f = None
            elif geo["lcal"] == "known_enter":
                f = None
            elif f is not None and geo["lcal"] == "known_pluck":
                T = geo["lcal_T"]
                L_eff = float(np.sqrt(T * G / (4.0 * mu * f * f)))
                propose_length(L_eff, [f"Known tension  {T:.0f} kgf", f"Measured note  {f:.0f} Hz", "",
                                       f"Working length = sqrt(T / (4 mu f^2)) = {L_eff*1000:.0f} mm"])
                set_overlay(cal_instructions())
            elif f is not None and geo["lcal"] == "fret_normal":
                geo["lcal_f"] = [f]; geo["lcal"] = "fret_clamped"
                det.fmin, det.fmax = fmin, fmax * geo["L"] / (args.fret_length / 1000.0) * 1.1
                set_overlay(cal_instructions())
            elif f is not None and geo["lcal"] == "fret_clamped":
                f_norm = geo["lcal_f"][0]
                det.fmin, det.fmax = band_normal
                if f <= f_norm * 1.03:
                    geo["lcal_msg"] = [f"Clamped note {f:.0f} Hz is not higher than the normal",
                                       f"note {f_norm:.0f} Hz. The crossing was not held still,",
                                       "or the pluck was on the hub side of the crossing.",
                                       "", "No change made.  Press l to try again."]
                    geo["lcal"] = "done"
                else:
                    L_eff = (args.fret_length / 1000.0) * f / f_norm
                    propose_length(L_eff, [f"Normal note   {f_norm:.0f} Hz", f"Clamped note  {f:.0f} Hz", "",
                                           f"Working length = {args.fret_length:.0f} x {f:.0f}/{f_norm:.0f} = {L_eff*1000:.0f} mm"])
                set_overlay(cal_instructions())
            elif f is not None and geo["lcal"] == "free":
                geo["lcal_f"] = [f]; geo["lcal"] = "crossed"
                det.fmin, det.fmax = band_normal
                set_overlay(cal_instructions())
            elif f is not None and geo["lcal"] == "crossed":
                f_free = geo["lcal_f"][0]
                if f <= f_free * 1.05:
                    geo["lcal_msg"] = [f"Crossed note {f:.0f} Hz is not higher than the",
                                       f"full-length note {f_free:.0f} Hz, so the crossing was",
                                       "probably still touching in step 1.",
                                       "", "No change made.  Press l to try again."]
                    geo["lcal"] = "done"
                else:
                    L_eff = geo["Lfull"] * f_free / f
                    propose_length(L_eff, [f"Full-length note  {f_free:.0f} Hz", f"Crossed note      {f:.0f} Hz", "",
                                           f"Working length = {geo['Lfull']*1000:.0f} x {f_free:.0f}/{f:.0f} = {L_eff*1000:.0f} mm"])
                set_overlay(cal_instructions())
            f = None
        if f is not None and st["armed"]:
            kgf = freq_to_kgf(f, geo["L"], mu)
            wheel.record(kgf, f, getattr(det, "accept_view", None))
            st["last"] = (wheel.cur + 1, f, kgf); st["last_t"] = now
            st["msg"] = f"spoke {wheel.cur+1}: {f:.0f} Hz = {kgf:.0f} kgf"
            if m_cur and n_cur >= 3 and abs(kgf - m_cur) / m_cur > 0.30:
                st["msg"] += f"   SUSPECT: {abs(kgf-m_cur)/m_cur:.0%} from side mean {m_cur:.0f} - re-pluck to confirm"
                wheel.suspect.add(wheel.cur)
            else:
                wheel.suspect.discard(wheel.cur)
            print(f"spoke {wheel.cur+1:2d} ({wheel.side_name(wheel.cur):8s})  {f:6.1f} Hz  {kgf:5.1f} kgf")
            if args.record:
                st["nrec"] += 1
                write_wav(f"{args.record}/spoke{wheel.cur+1:02d}_{st['nrec']:03d}.wav",
                          buf[-int(fs * 0.8):], fs)
            if wheel.calib is None:
                wheel.cur = (wheel.cur + 1) % wheel.n

        # ---- wheel map
        for i, (b_, v, lab) in enumerate(zip(bars, wheel.t, bar_labels)):
            b_.set_height(v if v is not None else 0)
            b_.set_hatch("//" if i in wheel.suspect else None)
            lab.set_y((v if v is not None else 0) + 1)
            lab.set_text((f"{v:.0f}?" if i in wheel.suspect else f"{v:.0f}") if v is not None else "")
        mL, sL, nL = wheel.side_stats("L"); mR, sR, nR = wheel.side_stats("R")
        mean_L.set_visible(mL is not None); mean_R.set_visible(mR is not None)
        if mL is not None: mean_L.set_ydata([mL, mL])
        if mR is not None: mean_R.set_ydata([mR, mR])
        cur_marker.set_x(wheel.cur + 0.5); cur_marker.set_width(1.0)

        # ---- stored spectrum of the selected spoke
        sp = wheel.spec[wheel.cur]
        if sp is not None:
            spec_stored.set_data(np.arange(len(sp)), sp)
            fc = wheel.f[wheel.cur]
            if fc:
                spec_fline.set_xdata([fc, fc]); spec_fline.set_visible(True)
                spec_label.set_text(f"spoke {wheel.cur+1}: {fc:.0f} Hz = {wheel.t[wheel.cur]:.0f} kgf")
        else:
            spec_stored.set_data([], []); spec_fline.set_visible(False)
            spec_label.set_text(f"spoke {wheel.cur+1}: no reading")

        # ---- reading panel
        if st["last"]:
            n_, f_, k_ = st["last"]
            read_big.set_text(f"{k_:.0f} kgf"); read_sub.set_text(f"spoke {n_}   {f_:.0f} Hz")
            read_big.set_color("black" if now - st["last_t"] < 4 else "0.55")
        else:
            read_big.set_text("--"); read_sub.set_text("")
        if geo["lcal"] is not None:
            state_txt.set_text("CALIBRATING LENGTH  (spoke readings paused)"); state_txt.set_color("tab:red")
        else:
            state_txt.set_text(("LISTENING" if st["armed"] else "PAUSED  (space to start)")
                               + f"     length {geo['L']*1000:.0f} mm"
                               + ("  UNCALIBRATED" if wheel.length_m is None else ""))
            state_txt.set_color("tab:green" if st["armed"] else "tab:red")

        # ---- learn turn gain from the last pass, once every spoke has been re-read
        if wheel.pass_log and mL is not None and mR is not None and nL + nR == wheel.n:
            pmL, pmR, tl, tr = wheel.pass_log[-1]
            if len(wheel.pass_log[-1]) == 4 and (tl or tr):
                g = []
                if tl: g.append((mL - pmL) / tl)
                if tr: g.append((mR - pmR) / tr)
                g = [v for v in g if v > 0]
                if g:
                    wheel.gain_samples.append(float(np.mean(g)))
                    wheel.gain = float(np.mean(wheel.gain_samples[-4:]))
                    wheel.pass_log[-1].append("used"); wheel.save()
                    st["msg"] = f"measured gain from last pass: {np.mean(g):.0f} kgf/turn (now using {wheel.gain:.0f})"

        # ---- raise-to-target plan
        plan_lines = []
        if args.target and mL is not None and mR is not None:
            tgt_L, tgt_R = args.target, args.target * wheel.ratio
            need_L, need_R = (tgt_L - mL) / wheel.gain, (tgt_R - mR) / wheel.gain
            tl = float(np.clip(need_L, -args.max_pass, args.max_pass))
            tr = float(np.clip(need_R, -args.max_pass, args.max_pass))
            # keep the sides moving in proportion so the dish doesn't walk
            if abs(need_L) > 1e-6 and abs(tl) < abs(need_L):
                tr = float(np.clip(need_R * tl / need_L, -args.max_pass, args.max_pass))
            tl, tr = round(tl * 8) / 8, round(tr * 8) / 8
            st["plan"] = (tl, tr)
            plan_lines = ["", f"RAISE TO TARGET  disc {tgt_L:.0f}  non-disc {tgt_R:.0f}",
                          f"  disc side needs {tgt_L-mL:+.0f} kgf (~{fmt_turns(need_L)} turn total)",
                          f"  non-disc needs  {tgt_R-mR:+.0f} kgf (~{fmt_turns(need_R)} turn total)"]
            if abs(tl) >= 1 / 8 or abs(tr) >= 1 / 8:
                plan_lines += [f"  THIS PASS: every disc nipple {fmt_turns(tl)} turn,",
                               f"             every non-disc nipple {fmt_turns(tr)} turn,",
                               f"  then press p and re-read all spokes"]
            else:
                plan_lines += ["  at target - equalise and true"]

        # ---- info panel
        lines = [f"current spoke: {wheel.cur+1}  ({wheel.side_name(wheel.cur)})"
                 + (f"   typing: {st['digits']}" if st["digits"] else ""), ""]
        if mL is not None:
            lines.append(f"disc      mean {mL:5.1f}  sd {sL:4.1f}  [{nL}/{wheel.n//2}]")
        if mR is not None:
            lines.append(f"non-disc  mean {mR:5.1f}  sd {sR:4.1f}  [{nR}/{wheel.n//2}]")
        if mL is not None and mR is not None:
            lines.append(f"ratio R/L {mR/mL:.2f}   (geometry {wheel.ratio:.2f})")
        lines.append(f"length {geo['L']*1000:.0f} mm" + (" (DEFAULT - press l!)" if wheel.length_m is None else " (cal)")
                     + f"   gain {wheel.gain:.0f} kgf/turn" + ("" if wheel.gain_samples else " (default)"))
        lines += plan_lines
        rec = [r for r in wheel.recommendations() if abs(r[3]) >= 1 / 16]
        if rec:
            lines += ["", "to equalise (biggest first):"]
            for i, v, d, turns in rec[:6]:
                lines.append(f" {i+1:2d} {wheel.side_name(i)[:4]:4s} {v:5.1f}  {fmt_turns(turns):>7} turn")
            if len(rec) > 6:
                lines.append(f" ... {len(rec)-6} more")
        elif nL + nR == wheel.n:
            lines += ["", "all spokes within 1/16 turn of side mean"]
        info_t.set_text("\n".join(lines))

        status.set_text(f"{st['msg']}    |    " + (f"detector: {det.last_reason}  (prom {det.last_prom:.1f})"
                                                 if st["armed"] else "mic off"))
        return ()

    ani = FuncAnimation(fig, update, interval=40, blit=False, cache_frame_data=False)
    try:
        plt.show()
    except KeyboardInterrupt:
        pass
    finally:
        shutdown()


if __name__ == "__main__":
    main()
