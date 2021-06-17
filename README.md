# nf-mapping

**My little DNA read mapping workflow which I wrote to use on sequences from museum specimen**

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04)
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10)
3. Download the pipeline, create a nextflow config profile that matches your cluster set-up ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) ) and start running your own analysis! If you want to use the existing `rackham.config` remember to specify your SNIC project ID (format: `snic20XX-XX-XXX`) as well as the path to `nf-polish/environment.yml`

    ```bash
    nextflow run nf-mapping/main.nf -profile rackham --reads 'READS'
    ```
4. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```
**Make sure** that there is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored).


