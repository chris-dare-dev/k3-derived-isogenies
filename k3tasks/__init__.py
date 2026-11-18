"""The computations behind the notes as tasks for the batch scheduler (taskbatch/).

`registry` lists the computations, `workspace` runs one away from the stored results, `compare` tells
whether a rerun reproduced the stored files, `tasks` connects them to the runner protocol, and `report`
adds the reproduction table to the scheduler's report.
"""

__version__ = "1.0"
