# spoke_tension

An acoustic spoke tension meter and wheel equaliser. Pluck a spoke near a USB microphone and the app:

1. finds the note,
2. converts it to tension (kgf),
3. stores it against a spoke index on a wheel map, and
4. recommends nipple turns to even out the wheel and bring it up to a target tension.

It is a single-file Python/Matplotlib desktop tool (`spoke_tension.py`).

Licence: GPL-3.0-or-later (see `LICENSE`).

## Install

Requirements: Python 3.9+ (developed on 3.14), a USB microphone, and a desktop session (the UI is a Matplotlib window).

### 1. System packages (Linux, Debian/Ubuntu)

```bash
sudo apt install python3-venv python3-tk libportaudio2
```

- `libportaudio2` is needed by `sounddevice` for audio input.
- `python3-tk` provides the Tk backend Matplotlib uses for the window. Other backends (Qt) work if you have them.
- macOS: `brew install portaudio`. Windows: nothing extra, the `sounddevice` wheel bundles PortAudio.

### 2. Clone and create the virtual environment

```bash
git clone git@github.com:sherwoac/spoke_tensioner.git
cd spoke_tensioner
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Python dependencies (`requirements.txt`): `numpy`, `sounddevice`, `matplotlib`. Last tested with numpy 2.5, sounddevice 0.5, matplotlib 3.11.

The examples below assume the environment is active. If you keep the venv elsewhere (for example `~/spoke_tensioner`), activate that one instead.

### 3. Find your microphone

```bash
python -c "import sounddevice as sd; print(sd.query_devices())"
```

Pass the index or a name substring as `--device` (for example `--device 0` or `--device rode`). If you get no input, check the OS lets the terminal use the microphone.

## Run

```bash
python spoke_tension.py --device rode --length 261 --full-length 271 --fret-length 188 \
    --target 115 --session ~/wheel_front.json
```

Analyse a single recording instead of using a live mic:

```bash
python spoke_tension.py --file pluck.wav
```

Useful options (see `--help` for all):

| Option | Meaning |
|---|---|
| `--device` | Input device (index or name substring) |
| `--length` | Working (acoustic) length in mm, see caveats below |
| `--full-length`, `--fret-length` | Whole free spoke / crossing-to-nipple length, used by calibration and detector guards |
| `--diameter` | Thin-section spoke diameter, default 1.8 mm |
| `--spokes` | Spoke count, default 28 |
| `--flange-left`, `--flange-right` | Hub flange distances in mm, used for the expected side ratio |
| `--target` | Disc-side target tension in kgf; non-disc target is scaled by the side ratio |
| `--max-pass` | Largest uniform adjustment per pass, default 1/4 turn |
| `--session` | JSON session file, auto-saved after every reading |
| `--record` | Save each accepted pluck as a wav |

## Keys

| Key | Action |
|---|---|
| `space` | Start/stop listening (the app starts paused) |
| `enter` / `→` | Next spoke |
| `←` | Previous spoke |
| digits then `enter` | Jump to spoke N |
| `u` | Undo last reading |
| `s` | Save session now |
| `c` | Calibrate turn gain: press, turn the nipple +1/4 turn, press again |
| `l` | Calibrate working length (three methods, each confirmed with y/n) |
| `p` | Log a uniform pass so gain is learned from the next read |
| `q` / `esc` / close window / Ctrl-C | Quit cleanly, session saved |

## Workflow

1. **Read** all spokes (the app auto-advances after each pluck).
2. **Equalise**: per-side means, and a list of spokes with the turns needed, in eighths.
3. **Raise**: with `--target`, plan uniform passes capped at `--max-pass`, keeping both sides proportional.

Spoke 1 is the spoke beside the valve on the disc side, and sides alternate around the rim.

## Example: front wheel

The wheel used to develop the tool: 28 spokes, 3x interlaced, Sapim Race 280 mm double-butted spokes (1.8 mm thin section), Chris King-style hub with the disc on the left.

### The lengths

A spoke has several lengths and only one of them is the right `L` for the formula:

```
 hub flange                    crossing                         rim
 (spoke head)                  with partner spoke               (nipple)
     |                              |                              |
     o==============================X==============================o
     |<----------------------------- 280 mm spoke ---------------->|
     |<---------------------- 271 mm free length ----------------->|   --full-length 271
     |                              |<------ 188 mm ------------->|   --fret-length 188
     |<---------------------- 261 mm working length ------------->|   --length 261
```

| Length | Value | Option | What it is |
|---|---|---|---|
| Spoke | 280 mm | | Manufacturer length, tip to elbow |
| Free length | 271 mm | `--full-length` | Whole free spoke, elbow to nipple seat, with the crossing partner held clear |
| Crossing to nipple | 188 mm | `--fret-length` | Distance from the crossing to the nipple, used by calibration method 2 |
| **Working length** | **261 mm** | `--length` | **The `L` in `T = 4 L² f² μ`.** Found by matching a phone app at about 80 kgf, not measured with a tape |

On this wheel the crossing is a soft node: it does not clamp the spoke, so the normal pluck tone behaves as if almost the whole spoke vibrates. Using the 188 mm tape measurement gives tensions about 40% too low. The working length of 261 mm is not yet verified against a tensiometer.

Hub geometry for the same wheel: flange distances 22.2 mm (left, disc) and 31.3 mm (right, non-disc), which are the defaults for `--flange-left` and `--flange-right`.

### The command

```bash
python spoke_tension.py --device 0 --length 261 --full-length 271 --fret-length 188 \
    --target 115 --session wheel_front.json
```

- `--device 0` is the USB microphone (index from `sd.query_devices()`).
- `--target 115` is the final tension for the disc side in kgf. The non-disc target is 115 × the side ratio.
- `--session wheel_front.json` stores readings and is reused next time you run, so you can stop and resume.

### The session file

[`examples/wheel_front.json`](examples/wheel_front.json) is a complete read of this wheel. The stored spectra are removed to keep it small (the real file is about 600 kB). The app loads it fine without them.

| Field | Meaning |
|---|---|
| `n` | Number of spokes (28) |
| `sides` | Side of each spoke by index: `L` = disc, `R` = non-disc. Spoke 1 is beside the valve and sides alternate |
| `t` | Tension in kgf for each spoke, in order |
| `f` | Measured pluck frequency in Hz for each spoke |
| `spec` | 1 Hz spectrum of each pluck (nulled in the example) |
| `length_m` | Working length used, in metres (`0.261`) |
| `ratio` | Expected non-disc/disc ratio from the flanges, 22.2 / 31.3 = 0.709 |
| `gain` | Kgf per full nipple turn used for turn recommendations, default 40 |
| `gain_samples` | Gain measurements learned from logged passes |
| `pass_log` | Uniform passes you logged with `p` |

To try it without a microphone, copy it and point the app at it:

```bash
cp examples/wheel_front.json my_wheel.json
python spoke_tension.py --session my_wheel.json --length 261 --target 115
```

### What the result looked like

| Side | Mean | Range | Std dev |
|---|---|---|---|
| Disc (L) | 115.0 kgf | 110.5 – 119.9 | 2.7 |
| Non-disc (R) | 94.5 kgf | 90.5 – 98.3 | 2.5 |

The measured ratio is 94.5 / 115.0 = 0.82, against 0.71 expected from the flanges. That may mean the rim sits about 2 mm toward the non-disc side, so check dish before trusting the ratio.

## Physics and caveats — please read

Tension is `T = 4 L² f² μ`, where `μ` is the mass per metre of the thin section of a double-butted spoke (about 20.1 g/m for 1.8 mm steel). In kgf this is roughly `0.0082 · L² · f²` with `L` in metres.

- **`L` is the working acoustic length, not a tape measurement.** On the wheel this was developed on (28H, 3x, interlaced, Sapim Race 280 mm) the spoke crossing acts as a soft node. The pluck tone behaves as if `L ≈ 261 mm`, although crossing-to-nipple is 188 mm and the free spoke is 271 mm. Using 188 under-reads by about 40%.
- **The 261 mm figure has not been verified against a tensiometer.** It came from agreement with a phone app at about 80 kgf. A longer `L` displays higher tension, so 261 is the safe-side choice. Calibrate against a known-tension reading (`l`, method 3) before trusting absolute numbers.
- **Treat the default length (186 mm) as wrong for most wheels.** The app warns when it is in use. Length priority is `--length`, then the session file, then `~/.spoke_tension_config.json`, then the 186 mm default.
- **Side ratio** (non-disc / disc) is expected near `flange_left / flange_right`. A measured ratio well off that may mean the rim is off-centre, so check dish first.

## Detection notes

Detection is deliberately conservative: it triggers on an abrupt onset, tracks the whole pluck until the tone fades, picks the pitch from a harmonic-weighted spectral peak in the tension band, and guards against octave errors. Readings more than 30% from the side mean are flagged suspect (hatched bar, "?") rather than rejected. Pausing aborts the audio stream entirely.

## Session files

Sessions are JSON (tensions, frequencies, spectra, gain samples, pass log, length). Saves are atomic: write a temp file, read it back and verify, rotate a `.bak`, then rename.

## Status and roadmap

Work in progress. Next up: known-tension calibration against a Park TM-1, dish and truing support, an offline test harness over labelled wav recordings, better quiet-pluck sensitivity, splitting into modules with tests, and packaging. The full list is in `CLAUDE.md`.
