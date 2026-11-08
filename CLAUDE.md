# spoke_tension — acoustic spoke tension meter and wheel equaliser

Single-file Python/Matplotlib desktop tool. Pluck a spoke near a USB mic; the app finds the
note, converts it to tension, stores it against a spoke index on a wheel map, and recommends
nipple turns to equalise and raise the wheel to a target. Licence: GPL-3.0-or-later.

Run (the venv lives at ~/spoke_tensioner):
    source ~/spoke_tensioner/bin/activate
    python spoke_tension.py --device rode --length 261 --full-length 271 --fret-length 188 \
        --target 115 --session ~/wheel_front.json

## Physics (don't change without re-reading this)
- T = 4 L² f² μ, μ = mass/metre of the thin (1.8 mm) section of a double-butted spoke
  (steel 7900 kg/m³ → 20.1 g/m). In kgf: T ≈ 0.0082 · L² · f² with L in metres.
- L is the **working (acoustic) length**, not a tape measurement. On this wheel (28H, 3x,
  interlaced, Sapim Race 280 mm) the crossing is a soft node: the normal pluck tone behaves as
  if L ≈ 261 mm even though crossing-to-nipple is 188 mm and the whole free spoke is 271 mm.
  Converting with 188 under-reads by ~40%. The 261 figure came from agreement with a
  commercial phone app at ~80 kgf; it has NOT been verified against a tensiometer.
  Error direction: a longer L displays HIGHER tension, so 261 is the safe-side choice.
- Length lives in three places, priority: `--length` > session file > ~/.spoke_tension_config.json
  > 186 mm default. Default is wrong for this wheel; the app warns loudly when it is in use.
- Side ratio: non-disc/disc expected ≈ flange_left/flange_right = 22.2/31.3 = 0.71
  (small-angle). Measured 0.82 at the end — may indicate rim ~2 mm off-centre toward the
  non-disc side; check dish before trusting the ratio.

## Detection pipeline (ToneDetector) — history of what broke and why it is this way
1. Energy-threshold triggers were abandoned: room noise + fixed threshold → false triggers
   and a stuck disarm. Now: abrupt onset (≥4× in-band energy AND rising edge AND a >6× step
   between consecutive 10 ms sub-frames). Speech ramps; typing has no sustained tone.
2. Judging the first 80 ms gave the attack transient, not the note. Now a pluck is TRACKED
   from onset until the tone fades (3 empty windows) or 0.7 s; one reading per pluck.
3. Per-window pitch: harmonic-weighted spectral peak (f + 0.7·|2f| + 0.4·|3f|) restricted to
   the tension band; autocorrelation on a narrow band around that peak only CONFIRMS a tone.
   Autocorrelation as a competitor was unreliable with two spokes ringing.
4. Guards: octave-up (shortest strong lag in AC), octave-down (strong peak near f/2 wins),
   full-length-mode (peak at L_full/L ratio → take the higher) — the last one is disabled when
   L_full/L < 1.15 because it misbehaves when the two lengths are similar.
5. Cluster choice across the track: score = max_amp · √windows; clusters >20% from the wheel
   prior (side mean, once ≥3 spokes read) are scaled ×0.3. Readings >30% from side mean are
   flagged SUSPECT (hatched bar, "?"), not rejected.
6. Pause (space) aborts the audio stream entirely; the app starts paused.

## Calibration (key l) — three methods, all bounds-checked and confirmed with y/n
1 free: hold crossing partner clear, pluck (L_full), then pluck normally → L = L_full·f_free/f.
   Only works at low tension. Gave ~240 mm early on.
2 fret: pluck normally, then pinch the crossing and brace against the rim so it cannot move,
   pluck → L = L_phys·f_clamped/f_normal. Pressing the two spokes together is NOT a node
   (note went DOWN); must be braced to the wheel. Unreliable so far.
3 known tension: type a tensiometer reading for the current spoke, pluck → L from physics.
   Most reliable. Not yet done on this wheel — do this before trusting absolute numbers.
Results outside [0.95·fret_length, 1.05·full_length] are rejected. A bad free-method
pluck once produced 161 mm and silently wrecked a session; that is why these guards exist.

## Session / persistence
- JSON session: tensions, frequencies, 1 Hz spectra per spoke, gain samples, pass log, length.
- Save is write-tmp → read back and verify → rotate .bak → rename. An interrupted non-atomic
  save once truncated the file to 0 bytes and lost a full wheel. Never regress this.
- Shutdown is one idempotent path (q / esc / window close / SIGINT) using stream.abort()
  (stop() blocks); a 200 ms GUI timer lets SIGINT through.

## Workflow the app supports
- Read all 28 (auto-advance; enter/arrows/N-enter to navigate; u undo).
- Equalise: per-side mean, list of spokes with turns in eighths (gain default 40 kgf/turn;
  c = quarter-turn calibration; p logs a uniform pass so gain is learned from the next read).
- Raise: --target sets disc target; non-disc target = target × ratio; plan caps uniform
  adjustment per pass at --max-pass (1/4 turn), keeps sides proportional.
- Hub: Chris King-style preload clamp on non-disc side; spoke 1 is beside the valve on the
  disc side; sides alternate around the rim.

## Open backlog (next things to do)
- [ ] Known-tension calibration against a Park TM-1 to pin L absolutely.
- [ ] Dish check / truing support: accept "rim pulls left at spokes a–b" and recommend paired
      adjustments; later maybe a dial indicator or camera.
- [ ] Offline test harness: run the detector over a folder of recorded .wav plucks labelled
      good/bad and report hits/misses/false alarms, so parameter changes are measured.
- [ ] Quiet-pluck sensitivity in noisy rooms (currently the main remaining miss case).
- [ ] Split the single file into modules (physics, detector, session, gui) and add tests.
- [ ] Package: pyproject.toml, console entry point, README with the physics caveats above.
- [ ] Linux: document `sudo apt install libportaudio2` for sounddevice.
