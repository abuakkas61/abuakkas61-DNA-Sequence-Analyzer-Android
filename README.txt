DNA SEQUENCE ANALYZER — ANDROID PROJECT (V12 ENGINE)

Contents
- main.py: Kivy Android interface
- sequence_engine.py: extracted V12 sequence and ORF engines
- buildozer.spec: Buildozer configuration

Features included
- IUPAC DNA/RNA sequence cleaning and validation
- Composition, GC/AT%, skew, purines/pyrimidines, entropy
- Reverse, complement, reverse complement, transcription
- Six-frame translation, codon counts
- ORF detection with adjustable minimum amino-acid length
- Sliding-window GC, GC/AT skew and entropy plots
- Positional sequence comparison (not alignment)
- FASTA open/save and text report export

Important
1. The original Tkinter file is not modified.
2. The sequence engine methods/constants were extracted from the supplied V12 source.
3. The Android interface is a Kivy rewrite. It does not use Tkinter or TkAgg.
4. This project has not been compiled into an APK in this environment. Test the app before building.
5. The file chooser initially browses the app's private data directory. Android shared-storage access differs by Android version; a Storage Access Framework integration may be needed for broad access to Downloads/Documents.

Build on Ubuntu / WSL2 (recommended)
  sudo apt update
  sudo apt install -y python3 python3-pip python3-venv git zip unzip openjdk-17-jdk       build-essential autoconf libtool pkg-config zlib1g-dev libncurses5-dev       libncursesw5-dev libtinfo6 cmake libffi-dev libssl-dev

  python3 -m venv ~/venvs/dna-apk
  source ~/venvs/dna-apk/bin/activate
  pip install --upgrade pip setuptools wheel
  pip install "Cython<3" buildozer

Copy/extract this project into WSL's Linux filesystem (for example ~/DNA_Sequence_Analyzer_Android), then:
  cd ~/DNA_Sequence_Analyzer_Android
  buildozer android debug

First build downloads Android SDK/NDK and python-for-android toolchain components and can take a while.
APK output should appear in bin/.

Desktop smoke test
  pip install kivy
  python main.py

If Buildozer fails, share the final 80-120 lines of the build log. Do not install python-for-android with pip as a workaround; Buildozer manages its toolchain.
