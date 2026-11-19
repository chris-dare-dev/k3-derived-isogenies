"""The computations behind the notes as tasks for the batch scheduler (taskbatch/).

`registry` lists the computations, `workspace` runs one away from the stored results, `compare` tells
whether a file of a rerun is the stored file, `expectation` holds the words of a result (what is
predicted, what was found, the verdict), `tasks` turns a batch into tasks and a task into its result for
the runner protocol, and `report` adds the reproduction table to the scheduler's report.
"""

__version__ = "1.0"
