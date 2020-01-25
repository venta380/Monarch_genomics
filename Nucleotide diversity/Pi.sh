#! /bin/bash -l
#SBATCH -A snic2019-3-35
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 30:00:00
#SBATCH -J Call_pi_all_sites
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL


python Monarch_pi_allsites.py 0 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 1 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 2 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 3 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 4 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
wait
python Monarch_pi_allsites.py 5 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 6 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 7 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 8 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
sleep 10m
python Monarch_pi_allsites.py 9 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi West &
wait



python Monarch_pi_allsites.py 0 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 1 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 2 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 3 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 4 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
wait
python Monarch_pi_allsites.py 5 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 6 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 7 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 8 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
sleep 10m
python Monarch_pi_allsites.py 9 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi East &
wait
