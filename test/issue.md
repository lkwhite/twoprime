I find the behavior of config file specification a bit unexpected.

If a `Snakefile` contains a `configfile` directive, it appears that this overrides a file specified on the command line with `--configfile`:

```
$ snakemake -npr -s pipeline/Snakefile --configfile test/config.yml 
WorkflowError in line 10 of /Users/jayhesselberth/devel/twoprime/pipeline/Snakefile:
Config file test/config.yml not found.
  File "/Users/jayhesselberth/devel/twoprime/pipeline/Snakefile", line 10, in <module>
```

The `Snakefile` has an identical `configfile` path, but `snakemake` is looking for the config file using  a `workdir` path that is also specified in the `Snakefile`. Specifying `--configfile` on the command line skips the `workdir` path requirement.

```
from datetime import date
today = date.today().isoformat()
workdir: 'results/%s' % today
```

However, when I try to tell `snakemake` where to find the actual config file, it throws the error above.

Does it seem reasonable to change this behavior to respect the specified `--configfile` on the command line over the one specified in the `Snakefile`? I can submit a PR if so.







