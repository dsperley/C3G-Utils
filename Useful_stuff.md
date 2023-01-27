---
editor_options: 
  markdown: 
    wrap: 72
---

[HPC](#HPC)\
[Bash](#Bash)\
[Git](#Git)\
[SSH](#SSH)\
[mysql](#mysql)\
[ftp](#ftp)\
[R](#r) [Other](#other)\
[Mysql](#mysql)

# HPC <a name="HPC"></a> {#hpc}

    $R_TOOLS = $MUGQIC_INSTALL_HOME/software/mugqic_tools/mugqic_tools-2.10.10/R-tools 

: z \## Jobs on abacus

    qsub -I -l walltime=2:00:00

## Scp from abacus

The key is -P 2222

    scp -r -P 2222 dperley@127.0.0.1:/lb/project/C3G/projects/Susta_Chicken_Breast_Myopathy_RNA-seq/rproj_Chicken_RNAseq .

## Interactive sessions

    salloc --time=2:00:0 --account=rrg-bourqueg-ad --mem=64000M
    salloc --time=6:0:0 --ntasks=40 --mem-per-cpu 2G  --account=def-bourqueg

## Location of R tools:

    /cvmfs/soft.mugqic/CentOS6/software/mugqic_tools/mugqic_tools-2.6.0/R-tools

## sbatch header

    #!/bin/bash
    #SBATCH --job-name=meth
    #SBATCH --account=rrg-bourqueg-ad
    #SBATCH --time=7:00:00
    #SBATCH --mem=100G
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=danielle.perley@mcgill.ca
    #SBATCH -n 1
    #SBATCH --output=meth_out_%j.txt

Or

    sbatch --job-name=meth_analysis --cpus-per-task=8 --account=rrg-bourqueg-ad --time=7:00:00 --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_%j.txt ./Methylation_analysis_walt.R Methylation_analysis_no_amb_reads | grep [0-9] | cut -f4

using an array

    #SBATCH --job-name=cutadapt
    #SBATCH --account=rrg-bourqueg-ad
    #SBATCH --time=20:00
    #SBATCH --mem-per-cpu=400M
    #SBATCH --cpus-per-task=5
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=danielle.perley@mcgill.ca
    #SBATCH --array=1-15
    #SBATCH -o ./logs/%A_%a.out

    module unload mugqic/python/2.7.14   
    module load python/3.9.6 mugqic/cutadapt/2.10 fastqc/0.11.9


    PROJECT="/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401"
    SHEET="${PROJECT}/Adapters.txt"
    OUT="${PROJECT}/02_Trim"

    if [ ! -d "$OUT" ]; then mkdir "$OUT"; fi


    adapter=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SHEET} | cut -f 2)
    sample

## to switch login nodes

    [dsperley@beluga2 ~]$ ssh beluga5

## get job efficency stats

    [dsperley@beluga2 Methylation_analysis]$ seff 21446876
    Job ID: 21446876
    Cluster: beluga
    User/Group: dsperley/dsperley
    State: FAILED (exit code 1)
    Cores: 1
    CPU Utilized: 01:44:07
    CPU Efficiency: 95.78% of 01:48:42 core-walltime
    Job Wall-clock time: 01:48:42
    Memory Utilized: 6.31 GB
    Memory Efficiency: 1.26% of 500.00 GB

## Monitoring active job

    sstat  -j 27190660 -a --format=JobID,MaxVMSize,MaxRSS
           JobID  MaxVMSize     MaxRSS 
    ------------ ---------- ---------- 
    27190660.ex+                       
    27190660.in+    288700K      7888K 

## Sacct

    [dsperley@cedar5 Zhenbao_Nanopore]$ sacct -j 28610945_2 --format=JobID%20,AllocCPUS,Elapsed,AllocTRES%40,MaxRSS,MaxVMSize,State --unit=G
                   JobID  AllocCPUS    Elapsed                                AllocTRES     MaxRSS  MaxVMSize      State 
    -------------------- ---------- ---------- ---------------------------------------- ---------- ---------- ---------- 
              28610945_2          8   20:15:09       billing=13,cpu=8,mem=56000M,node=1                       CANCELLED+ 
        28610945_2.batch          8   20:15:15                  cpu=8,mem=56000M,node=1   5593560K    155692K  CANCELLED 
       28610945_2.extern          8   20:15:15       billing=13,cpu=8,mem=56000M,node=1          0      4364K  COMPLETED 

## checking fair share

    [dsperley@cedar1 rnaseqRep]$ sshare -l -u dsperley | grep rrg-bourqueg-ad
     rrg-bourqueg-ad_cpu                162000    0.001215   227503923    0.002591      0.002591              0.469100                                cpu=904477,mem=3965224853,ene+ 
      rrg-bourqueg-ad_c+   dsperley          1    0.017544    20963900    0.000239      0.092148   0.357268   0.190389                                cpu=0,mem=0,energy=0,node=0,b+ 

## installing R packages on CC

Use module load r/4.1.0 instead module load
mugqic/R_Bioconductor/4.1.0_3.13

## cpus

from slack

@po I am did modification in many pipelines to have --mem-per-cpu. It
has many advantages. Also I would recommend to use M instead of G.
Beluga has 191000M and 95000M and 771000M machines. --mem-per-cpu value
that are right are: 4775M, 2375M, and 19275M. All other value will
potentially make your job wait longer on queue.

@Mathieu FYI: A rule of thumb for reads of \~100bp is to set
MAX_RECORDS_IN_RAM to be 250,000 reads per each GB given to the -Xmx
parameter (picard)

## Transferring between file systems

(<https://docs.computecanada.ca/wiki/Frequently_Asked_Questions#Moving_files_from_project_to_scratch_or_home_filesystems>)

Moving files from project to scratch or home filesystems If you want to
move files from your project into your scratch or home space, you should
not use the mv command. Instead, we recommend using the regular cp, or
the rsync command.

It is very important to run cp and rsync correctly to ensure that the
files copied over to the project space have the correct group ownership.
With cp, do not use the archive -a option. And when using rsync, make
sure you specify the --no-g --no-p options, like so:

    [name@server ~]$ rsync -axvH --no-g --no-p  $HOME/projects/<project>/some_other_directory $HOME/scratch/some_directory

## sbatch

    cat validate_sam.sh | sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID  -o validate.out -J validate.Sam --time=24:00:00 --mem 12G -c 1 -N 1 -q centos7  

## archiving

    dar -v -c /nearline/rrg-vmooser/dsperley/pad_exome_plus_fastqs -[ fastq_files.txt -s 300G

## cancel all jobs

    sq | awk '{print $1}' | tail -n +2 | xargs scancel 

## mounts

If you mount abacus through sshfs (in mac) and if it gives input/output
error when you sleep the mac and come back, use this workaround to
unmount the drives. diskutil umount force [drive path]; if this doesn't
work ps -ax \| grep sshfs to get the PID of sshfs and then kill the
sshfs process using PID kill 9 PID then umount -f [drive path]

[dperley\@Alyssum](mailto:dperley@Alyssum){.email} rrg-bourqueg-ad %
sshfs -o follow_symlinks -o IdentityFile=\~/.ssh/id_rsa -o
defer_permissions
[dsperley\@cedar.computecanada.ca](mailto:dsperley@cedar.computecanada.ca){.email}:
cedar [dperley\@Alyssum](mailto:dperley@Alyssum){.email} mnt % sshfs
[dsperley\@cedar.computecanada.ca](mailto:dsperley@cedar.computecanada.ca){.email}:
cedar

## list all accounts

sacctmgr show associations where user=dsperley
format=cluster,account%40,qos

# Bash <a name="Bash"></a> {#bash}

## transpose

    awk '/^CATEGORY/ {split($0,header);n=1;next; } {if(n!=1) next; for(i=2;i<=NF;++i) printf("%s\t%s\t%s\n",$1,header[i],$i);}' |\
    column -t

## get non - size 0

    find . ! -size 0 

## get the total disk usage from a list of files

    cat ~/BQC19_r3/BQC19/anonymisation/cram_filenames.txt | xargs -d \\n du -ch

## select 2 lines from file

    sed -n "2p;8p" GSP3_PAH_Q7_stats.txt 

## link specific fastq files

    for i in $(cut -f38 ../2022_16S_samples.txt); do  ln -s ../../2022/raw_fastq_dl/${i}_R1.fastq.gz .; echo ln -s ../../2022/raw_fastq_dl/${i}_R2.fastq.gz .; done

## find examples

    find . -type f -name "*sh" ! -path "*genpipe*" 
    find . -type f -name "*sh" ! -path "*genpipe*" -exec cp '{}' scripts_to_dl \;
    find . -maxdepth 1 -user "dsperley"
    find ~/projects/ctb-bourqueg-ac/HostSeqQc/GLOBUS_SHARE/BQC19/Complete_Clinical_And_Genomics/HostSeq/WGS -name "*cram" -printf "%f\t%CY-%Ch\n

    find raw_reads -type f -name "*I[12]_001.fastq.gz" -exec rm {} +

    find . -type f  ! -newermt 2022-08-18 -exec mv {} ../locatit_0817 \;
    find ./genpipes/job_output/bwa_mem_sambamba_sort_sam -name "*done" | grep -f files_to_redo_08152022.txt | xargs rm

    find genpipes/alignment -mindepth 1 -maxdepth 1 -type d -type d '!' -exec sh -c 'ls -1 "{}"|egrep -i -q ".*[0-9].properties$"' ';' -print 

    find . -maxdepth 1 -type d -newermt "2022-10-31" -name "BQC*" -exec mkdir -p {}/cram \;

    find /home/dsperley/projects/ctb-bourqueg-ac/HostSeqQc/GLOBUS_SHARE/BQC19/Complete_Clinical_And_Genomics/Release10_HS_Release11/WGS -type f | grep -f Hostseq_ids_to_remove.txt | xargs rm

## hard link multiple directories

    cp -Rl ../../Release7_APR/WGS/BQC* .

ls *cram* \| grep -f ../release7_batch1_legit_cram.tsv \| xargs rm

## awk examples

extract filename from path

    awk '{n=split($NF,a,"/"); print $1,a[n]}' WD_PD1_4_1_DES1351A292-38_S167_L001_R2_001.fastq.gz.md5

## rsync

    rsync --no-relative --no-dirs --files-from=plates4and5_files_to_copy.txt / dsperley@cedar.computecanada.ca:/project/6007512/C3G/projects/Desharnais_16S/plates4and5

## units

### ls

default for files (bytes), directory totals (kilobytes) du (kilobytes)

## vim

to remove lines containg a string

    :g/--knownAlleles/d
    :g/^\s*\\$/d
    \s = whitespace
    \\ backslach
    ## removes lines with only a slash

# Git <a name="Git"></a> {#git}

## checkout a remote branch

    git checkout -b importSpreadsheetData origin/importSpreadsheetData

## show remotes

    dperley@Alyssum megadata % git remote -v
    origin  https://dsperley1@bitbucket.org/genap/megadata.git (fetch)
    origin  https://dsperley1@bitbucket.org/genap/megadata.git (push)

## retrieve newer files on remote

    dperley@Alyssum megadata % git fetch origin
    remote: Enumerating objects: 103, done.
    remote: Counting objects: 100% (103/103), done.
    remote: Compressing objects: 100% (79/79), done.
    remote: Total 81 (delta 50), reused 2 (delta 0), pack-reused 0
    Unpacking objects: 100% (81/81), 697.29 KiB | 589.00 KiB/s, done.
    From https://bitbucket.org/genap/megadata
       faf1ea7..84c35cf  importSpreadsheetData -> origin/importSpreadsheetData

## Show information about remote

    git remote show origin
    git pull

# SSH <a name="SSH"><a/> {#ssh}

## how to add ssh-keys

(<https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent>)

Generating a new SSH key Open Terminal.

Paste the text below, substituting in your GitHub email address.

    $ ssh-keygen -t ed25519 -C "your_email@example.com"

Note: If you are using a legacy system that doesn't support the Ed25519
algorithm, use:

    $ ssh-keygen -t rsa -b 4096 -C "your_email@example.com"

This creates a new ssh key, using the provided email as a label.

    > Generating public/private ed25519 key pair.

When you're prompted to "Enter a file in which to save the key," press
Enter. This accepts the default file location.

    > Enter a file in which to save the key (/Users/you/.ssh/id_ed25519): [Press enter]

At the prompt, type a secure passphrase. For more information, see
"Working with SSH key passphrases."

    > Enter passphrase (empty for no passphrase): [Type a passphrase]
    > Enter same passphrase again: [Type passphrase again]

Adding your SSH key to the ssh-agent Before adding a new SSH key to the
ssh-agent to manage your keys, you should have checked for existing SSH
keys and generated a new SSH key. When adding your SSH key to the agent,
use the default macOS ssh-add command, and not an application installed
by macports, homebrew, or some other external source.

Start the ssh-agent in the background.

    $ eval "$(ssh-agent -s)"
    > Agent pid 59566

Depending on your environment, you may need to use a different command.
For example, you may need to use root access by running sudo -s -H
before starting the ssh-agent, or you may need to use exec ssh-agent
bash or exec ssh-agent zsh to run the ssh-agent.

If you're using macOS Sierra 10.12.2 or later, you will need to modify
your \~/.ssh/config file to automatically load keys into the ssh-agent
and store passphrases in your keychain.

First, check to see if your \~/.ssh/config file exists in the default
location.

    $ open ~/.ssh/config
    > The file /Users/you/.ssh/config does not exist.

If the file doesn't exist, create the file.

    $ touch ~/.ssh/config (permissions should be 644)

Open your \~/.ssh/config file, then modify the file to contain the
following lines. If your SSH key file has a different name or path than
the example code, modify the filename or path to match your current
setup.

    Host *
      AddKeysToAgent yes
      UseKeychain yes
      IdentityFile ~/.ssh/id_ed25519

Note: If you chose not to add a passphrase to your key, you should omit
the UseKeychain line.

Note: If you see an error like this

    /Users/USER/.ssh/config: line 16: 

Bad configuration option: usekeychain add an additional config line to
your Host \* section:

    Host *
      IgnoreUnknown UseKeychain

Add your SSH private key to the ssh-agent and store your passphrase in
the keychain. If you created your key with a different name, or if you
are adding an existing key that has a different name, replace id_ed25519
in the command with the name of your private key file.

    $ ssh-add -K ~/.ssh/id_ed25519

Note: The -K option is Apple's standard version of ssh-add, which stores
the passphrase in your keychain for you when you add an ssh key to the
ssh-agent. If you chose not to add a passphrase to your key, run the
command without the -K option.

If you don't have Apple's standard version installed, you may receive an
error. For more information on resolving this error, see "Error:
ssh-add: illegal option -- K."

Add the SSH key to your account on GitHub. For more information, see
"Adding a new SSH key to your GitHub account."

# mysql <a name="mysql"></a> {#mysql}

To start:

    brew services start mysql

log into my sql as root

     mysql -u root -p

View users
`mysql> select user,host from mysql.user; +------------------+-----------+ | user             | host      | +------------------+-----------+ | mysql.infoschema | localhost | | mysql.session    | localhost | | mysql.sys        | localhost | | root             | localhost | +------------------+-----------+ 4 rows in set (0.00 sec)`

\# FTP <a name="ftp"></a>

\- type passwords (no copying)

    sftp -c aes128-cbc vsoleimani@sftp.mcgillgenomecentre.ca 

remove directory

    rmdir dir

uploading all fastq files (do in a tmux session!)

    put *gz

downloading directories

    get -r *

Copying a directory via SFTP is slightly more complicated. This is
because most command line implementations of SFTP cannot directly copy a
directory, instead you can only copy the *contents* of a directory. In
practice, this means the outermost layer (i.e., the directory itself)
will not be copied. You will usually have to make it yourself.

To copy a directory, follow these steps:

1.  Create a local directory of your choice using the command *lmkdir
    directoryname*.

2.  Use *lcd directoryname* to navigate your working directory into the
    empty directory you\'ve just created.

3.  Move inside the remote directory you want to copy.

4.  From within the remote directory, copy all the files using the
    command *get -r \**.

This will copy all files and sub-directories contained within the
directory. This process works same way for transferring a directory to
the remote host, but uses the *put* command.

Another option is to simply compress the directory you\'re trying to
move, at which point the compressed folder can be transferred like any
other file.

# R <a name="R"></a> {#r}

## Installation/update

-   Download pkg from (not arm64)
    <https://cran.r-project.org/bin/macosx/>

install all packages from old version go to

    /Library/Frameworks/R.framework/Versions/4.0/Resources/library

where 4.0 is version and

    ls > ~/R_4.0.5_packages.txt

In R:

    old_packages <- read.delim("~/R_4.0.5_packages.txt")
    current <- installed.packages()
    to_install <- old_packages[!old_packages$V1 %in% current,]
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(version = "3.14")
    BiocManager::install(to_install)

    packages ‘EPIC’, ‘MCPcounter’, ‘SeuratData’, ‘SomaDataIO’, ‘hierarchicell’, ‘immunedeconv’, ‘pbmc3k.SeuratData’, ‘scPower’, ‘sfit’, ‘translations’, ‘xCell’ are not available for Bioconductor version '3.14'

# Other <a name="other"></a> {#other}

## installing genome

      ./install_genome.sh Delphinapterus_leucas Beluga ASM228892v3 None Ensembl 105

## Genap

    https://genap.ca/p/account/home

Can add another workspap for applications, need filebrower and datahub
to transfer files, in either datahub or file brower, add ssh key, found
in \~/.ssh/id_rsa.pub in home directory of system you want to transfer
from. ie cedar place keep in authorized keys file and save file now on
system thay you're transferring from:

pu9 is workspace id

    [dsperley@cedar1 Alignments]$ sftp -P 22004 sftp_pu9@sftp-arbutus.genap.ca
    Enter passphrase for key '/home/dsperley/.ssh/id_rsa': 
    Connected to sftp-arbutus.genap.ca.
    sftp> ls
    datahub641   share        
    sftp> cd datahub641/
    sftp> lls
    GSP3_Hind3_Q7.sam     GSP3_Hind3_Q7_sorted.bam.bai  GSP3_PAH_Q7.sam     GSP3_PAH_Q7_sorted.bam.bai
    GSP3_Hind3_Q7_sorted.bam  GSP3_Hind3_Q7_stats.txt   GSP3_PAH_Q7_sorted.bam  GSP3_PAH_Q7_stats.txt
    sftp> put *sorted*
    Uploading GSP3_Hind3_Q7_sorted.bam to /datahub641/GSP3_Hind3_Q7_sorted.bam
    GSP3_Hind3_Q7_sorted.bam                                                 100% 2200MB 101.8MB/s   00:21    
    Uploading GSP3_Hind3_Q7_sorted.bam.bai to /datahub641/GSP3_Hind3_Q7_sorted.bam.bai
    GSP3_Hind3_Q7_sorted.bam.bai                                             100% 1874KB  22.7MB/s   00:00    
    Uploading GSP3_PAH_Q7_sorted.bam to /datahub641/GSP3_PAH_Q7_sorted.bam
    GSP3_PAH_Q7_sorted.bam                                                   100% 2548MB 107.4MB/s   00:23    
    Uploading GSP3_PAH_Q7_sorted.bam.bai to /datahub641/GSP3_PAH_Q7_sorted.bam.bai
    GSP3_PAH_Q7_sorted.bam.bai   
     sftp> bye 

## GLOBUS command line

used instructions from ComputeCanada to install globus -cli
<https://docs.alliancecan.ca/wiki/Globus#Command_Line_Interface_.28CLI.29>

Create a virtual environment to install the Globus CLI into (see
creating and using a virtual environment). \## need a non mugqic python
(remove any mugqic python) \$ virtualenv \$HOME/.globus-cli-virtualenv
Activate the virtual environment \$ source
\$HOME/.globus-cli-virtualenv/bin/activate Install Globus CLI into the
virtual environment (see installing modules). \$ pip install globus-cli
Then deactivate the virtual environment. \$ deactivate To avoid having
to load that virtual environment every time before using Globus, you can
add it to your path. \$ export
PATH=$PATH:$HOME/.globus-cli-virtualenv/bin \$ echo 'export
PATH=$PATH:$HOME/.globus-cli-virtualenv/bin'\>\>\$HOME/.bashrc

then authenticated <https://docs.globus.org/cli/quickstart/>

end points 6c66d53d-a79d-11e8-96fa-0a6d4e044368 \|
[mcgilluniversity\@globusid.org](mailto:mcgilluniversity@globusid.org){.email}
\| mcgilluniversity#genomecentre-general
[[dsperley\@beluga3](mailto:dsperley@beluga3){.email} \~]\$ globus
endpoint search "computecanada" ID \| Owner \| Display Name\
------------------------------------ \|
------------------------------------------------------------ \|
--------------------------------- 278b9bfe-24da-11e9-9fa2-0a06afd4a22e
\|
[computecanada\@globusid.org](mailto:computecanada@globusid.org){.email}
\| computecanada#beluga-dtn\
c99fd40c-5545-11e7-beb6-22000b9a448b \|
[computecanada\@globusid.org](mailto:computecanada@globusid.org){.email}
\| computecanada#cedar-dtn\
dabdce63-6d04-11e5-ba46-22000b92c6ec \|
[computecanada\@globusid.org](mailto:computecanada@globusid.org){.email}
\| computecanada#mammouth\
77506016-4a51-11e8-8f88-0a6d4e044368 \|
[computecanada\@globusid.org](mailto:computecanada@globusid.org){.email}
\| computecanada#niagara\
4ccf5768-fcad-4cbb-a2c7-563a109c8bd3 \|
[4ccf5768-fcad-4cbb-a2c7-563a109c8bd3\@clients.auth.globus.org](mailto:4ccf5768-fcad-4cbb-a2c7-563a109c8bd3@clients.auth.globus.org){.email}
\| ComputeCanada - Narval

## Loading files from Genap to IGV

-   Key is for datahub to be set to public
-   then in IGV go File -\> load from url
-   copy url of each bam file

## submit genpipes

This will allow submit upto 500 jobs at once.

# Mysql <a name="mysql"></a>

## install

    brew install mysql
    "...."

    ==> mysql
    We've installed your MySQL database without a root password. To secure it run:
        mysql_secure_installation

    MySQL is configured to only allow connections from localhost by default

    To connect run:
        mysql -uroot

    To restart mysql after an upgrade:
      brew services restart mysql
    Or, if you don't want/need a background service you can just run:
      /opt/homebrew/opt/mysql/bin/mysqld_safe --datadir=/opt/homebrew/var/mysql

To start service

    dperley@Alyssum sql %  brew services start mysql
    ## password is Hobbes
    dperley@Alyssum sql % mysqladmin -u root password Hobbes
    ``
    create user

use mysql; CREATE USER 'dsperley'\@'localhost' IDENTIFIED BY 'Hobbs1!M';
GRANT ALL ON *.* TO 'dsperley'\@'localhost';

# Samtools and BCFtools

## get list of sample from vcf file

    bcftools query -l chr8.dose.vcf.gz

renaming chrm for i in {1..22} Y X M; do echo -e "$i\tchr$i"
\>\>chrNames.txt; done bcftools annotate --rename-chrs chrNames.txt
Homo_sapiens.GRCh38.dbSNP150.vcf.gz \| bgzip \>
Homo_sapiens.GRCh38.dbSNP150_renamed.vcf.gzv

awk -F "\t" '{print \$1,\$2-50,\$3+49}' snps.bed

## other bits and pieces

10-01-2012, 12:07 PM Thanks for responding. I finally found a solution
to this, so I'll post it here in case is useful for others.

I found out that samtools filters reads before including them in the
pileup; it reads the flag field in the bam file and discards reads that
a) are not paired b) not properly mapped c) mate is not mapped d)
alignment is not primary e) reads fail quality control of vendor f) is
marked as PCR duplicates.

If the filters (a) and (c) are not desired, you can use the parameter
-A. In addition, samtools performs realignment unless the parameter -B
is used, and discards low quality reads unless -Q0 is used. Finally, it
stops at a certain number of reads unless the -d parameter is invoqued.
