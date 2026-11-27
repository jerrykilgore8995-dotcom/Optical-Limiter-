"""Export the ngspice DC state of the unmodified manual circuit."""
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

source = Path(__file__).resolve().parents[1] / 'manual' / 'manual.cir'
dest = Path(sys.argv[1]).resolve()
dest.parent.mkdir(parents=True, exist_ok=True)
deck = re.sub(r'(?im)^\.include\s+(\S+)',
              lambda m: f'.include "{(source.parent / m[1]).resolve().as_posix()}"',
              source.read_text())
deck = re.sub(r'(?im)^\.(op|tran|end)\b[^\n]*', '', deck)
with TemporaryDirectory() as directory:
    work = Path(directory)
    (work / 'seed.cir').write_text(deck + '''
.control
set filetype=ascii
op
write seed.raw all
quit
.endc
.end
''')
    result = subprocess.run([os.environ.get('NGSPICE', 'ngspice'), '-b', 'seed.cir'],
                            cwd=work, capture_output=True, text=True, check=True)
    data = (work / 'seed.raw').read_text()
    if 'Plotname: Operating Point' not in data or 'No. Points: 1' not in data:
        raise RuntimeError(result.stdout + result.stderr)
    dest.write_text(data)
