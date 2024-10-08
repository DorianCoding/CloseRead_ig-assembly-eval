#!/bin/sh
# properties = {"type": "single", "rule": "coverageAnalysis", "local": false, "input": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/gene_position/test.final.Ig_loci.txt", "/home1/zhuyixin/zhuyixin_proj/AssmQuality/aligned_bam/test/test_merged_sorted_primary.bam", "/home1/zhuyixin/zhuyixin_proj/AssmQuality/aligned_bam/test/test_merged_sorted_primary.bam.csi", "code/coverage_snake.sh"], "output": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/errorStats/test/pileup.end"], "wildcards": {"HOME": "/home1/zhuyixin/zhuyixin_proj/AssmQuality", "species": "test"}, "params": {"species": "test", "assemblies": "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/test.merged.fasta", "IGH_out": "/home1/zhuyixin/zhuyixin_proj/AssmQuality/errorStats/test/IGH_pri_pileup.txt", "IGK_out": "/home1/zhuyixin/zhuyixin_proj/AssmQuality/errorStats/test/IGK_pri_pileup.txt", "IGL_out": "/home1/zhuyixin/zhuyixin_proj/AssmQuality/errorStats/test/IGL_pri_pileup.txt", "conda": "/spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh"}, "log": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/log/coverageAnalysis.test.log"], "threads": 10, "resources": {"mem_mb": 30000, "mem_mib": 46570, "disk_mb": 48832, "disk_mib": 46570, "tmpdir": "<TBD>", "mem": "30G"}, "jobid": 6, "cluster": {}}
cd /home1/zhuyixin/ig-assembly-eval && /home1/zhuyixin/.conda/envs/assembly/bin/python3.10 -m snakemake --snakefile '/home1/zhuyixin/ig-assembly-eval/Snakefile' --target-jobs 'coverageAnalysis:HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality,species=test' --allowed-rules 'coverageAnalysis' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=30000' 'mem_mib=46570' 'disk_mb=48832' 'disk_mib=46570' --wait-for-files '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.nxiqn1t0' '/home1/zhuyixin/zhuyixin_proj/AssmQuality/gene_position/test.final.Ig_loci.txt' '/home1/zhuyixin/zhuyixin_proj/AssmQuality/aligned_bam/test/test_merged_sorted_primary.bam' '/home1/zhuyixin/zhuyixin_proj/AssmQuality/aligned_bam/test/test_merged_sorted_primary.bam.csi' 'code/coverage_snake.sh' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'params' 'mtime' 'software-env' 'input' 'code' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 60000 --scheduler 'ilp' --scheduler-solver-path '/home1/zhuyixin/.conda/envs/assembly/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.nxiqn1t0/6.jobfinished' || (touch '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.nxiqn1t0/6.jobfailed'; exit 1)

