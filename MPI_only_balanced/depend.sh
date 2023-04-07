ID=($(sbatch send_job_compile.sh))
sbatch --dependency=afterok:${ID[3]} send_job_n1.sh
sbatch --dependency=afterok:${ID[3]} send_job_n2.sh
sbatch --dependency=afterok:${ID[3]} send_job_n3.sh
sbatch --dependency=afterok:${ID[3]} send_job_n4.sh
sbatch --dependency=afterok:${ID[3]} send_job_n6.sh
sbatch --dependency=afterok:${ID[3]} send_job_n8.sh
