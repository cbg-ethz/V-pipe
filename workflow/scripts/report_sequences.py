#!/usr/bin/env python

import hashlib
import json
import pathlib
import shutil
import sys
import tempfile
import zipfile
from datetime import datetime, timezone

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import crc64, seguid

LOCAL_TIMEZONE = datetime.now(timezone.utc).astimezone().tzinfo

results = []
raw_data_files = []

for f in sys.argv[1:]:
    f = pathlib.Path(f)
    _, sample, batch, *_ = f.parts

    timestamp = datetime.fromtimestamp(f.stat().st_ctime, tz=LOCAL_TIMEZONE)
    sequences = []

    dehuamized_raw_reads = f.parent.parent / "raw_data" / "dehuman.cram"
    if dehuamized_raw_reads.exists():
        raw_data_files.append(dehuamized_raw_reads)
        raw_data = str(dehuamized_raw_reads)
    else:
        raw_data = None

    for record in SeqIO.parse(f, "fasta"):
        sequence = record.seq
        sequences.append(
            dict(
                header=record.id,
                seguid=seguid(sequence),
                crc64=crc64(sequence),
                sequence=str(sequence),
            )
        )
    results.append(
        dict(
            sample=sample,
            batch=batch,
            created=str(timestamp),
            sequences=sequences,
            raw_data=raw_data,
        )
    )


folder = pathlib.Path(tempfile.mkdtemp())

with open(folder / "summary.json", "w") as fh:
    json.dump(results, fh, indent=4)

for p in raw_data_files:
    target = folder / p
    target.parent.mkdir(parents=True)
    print(p, target)
    shutil.copy(p, target)

with zipfile.ZipFile("summary.zip", "w", zipfile.ZIP_DEFLATED) as zip_file:
    for entry in folder.rglob("*"):
        zip_file.write(entry, entry.relative_to(folder))
