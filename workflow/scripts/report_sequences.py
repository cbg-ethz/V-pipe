#!/usr/bin/env python

import hashlib
import pathlib
import json
import sys
from datetime import datetime, timezone

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid, crc64

LOCAL_TIMEZONE = datetime.now(timezone.utc).astimezone().tzinfo

results = []

for f in sys.argv[1:]:
    f = pathlib.Path(f)
    _, sample, batch, *_ = f.parts

    timestamp = datetime.fromtimestamp(f.stat().st_ctime, tz=LOCAL_TIMEZONE)
    sequences = []
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
        dict(sample=sample, batch=batch, created=str(timestamp), sequences=sequences)
    )

with open("summary.json", "w") as fh:
    json.dump(results, fh, indent=4)
