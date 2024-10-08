#!/bin/sh
# properties = {"type": "single", "rule": "lociLocation", "local": false, "input": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/test.pri.fasta", "code/igDetective.sh"], "output": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/test.pri.txt", "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/test.alt.txt"], "wildcards": {"HOME": "/home1/zhuyixin/zhuyixin_proj/AssmQuality", "species": "test"}, "params": {"pri_outdir": "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/test.pri.igdetective/", "alt_outdir": "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/test.alt.igdetective/", "species": "test", "igdetective_home": "/home1/zhuyixin/IgDetective", "conda": "/spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh", "haploid": "True", "condaEnv": "/home1/zhuyixin/.conda/envs"}, "log": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/log/igDetective.test.log"], "threads": 10, "resources": {"mem_mb": 30000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "mem": "30G"}, "jobid": 5, "cluster": {}}
cd /home1/zhuyixin/ig-assembly-eval && /home1/zhuyixin/.conda/envs/assembly/bin/python3.10 -m snakemake --snakefile '/home1/zhuyixin/ig-assembly-eval/Snakefile' --target-jobs 'lociLocation:HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality,species=test' --allowed-rules 'lociLocation' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=30000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.b17l9uhn' '/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/test.pri.fasta' 'code/igDetective.sh' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'software-env' 'params' 'code' 'input' 'mtime' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 60000 --scheduler 'ilp' --scheduler-solver-path '/home1/zhuyixin/.conda/envs/assembly/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.b17l9uhn/5.jobfinished' || (touch '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.b17l9uhn/5.jobfailed'; exit 1)

