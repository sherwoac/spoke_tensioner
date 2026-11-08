# spoke_tension

An acoustic spoke tension meter and wheel equaliser. Pluck a spoke near a USB microphone and the app:

1. finds the note,
2. converts it to tension (kgf),
3. stores it against a spoke index on a wheel map, and
4. recommends nipple turns to even out the wheel and bring it up to a target tension.

It is a single-file Python/Matplotlib desktop tool (`spoke_tension.py`).

Licence: GPL-3.0-or-later (see `LICENSE`).

## Install

```bash
# Linux: PortAudio is needed by sounddevice
sudo apt install libportaudio2

python3 -m venv ~/spoke_tensioner
source ~/spoke_tensioner/bin/activate
pip install sounddevice numpy matplotlib
```

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
