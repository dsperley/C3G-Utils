!/bin/bash

echo "Source path= "$1

#/sb/send_offline/mugqic/projects
####################################
#
# Backup to NFS mount script.
#
####################################

# What to backup.
backup_foldr_n_path=$1
backup_foldr=`echo $backup_foldr_n_path | awk -F"/" '{print $NF}'`
echo
echo "Archive folder=" $backup_foldr
echo
echo "Archive directory path=" $backup_foldr_n_path

# Where to backup to.
dest_path="/sb/send_offline/mugqic/projects"

# Create archive filename.
date=$(date +%Y-%m-%d)
run_prefix="${backup_foldr}_${date}"
#tree $1 > _source_file_tree.txt
# Print start status message.
echo
echo "Backing up $1 to $dest_path/$backup_foldr"
date

rsync  -axvH --no-g --no-p --exclude="trim" --exclude="alignment*"  $1 $dest_path > ~/${run_prefix}.filtered.log
if [[ $? -gt 0 ]]
then
  # take failure action here
  echo "*******Rsync FAILED*************"
else
  echo "Rsync is Done"
# Backup the files using tar.
  cd  $dest_path
  tar -czvf ${backup_foldr}.tar.gz $backup_foldr > ~/${run_prefix}.tar.log
  echo -e '\nHit [Ctrl]+[D] to exit this child shell.'
  $SHELL
fi
rsync --dry-run -axvH --no-g --no-p  $1 $dest_path/$backup_foldr > ~/${run_prefix}.all.log
echo "Copied file  log in ~/${run_prefix}.filtered.log"
echo "All files og in ~/${run_prefix}.all.log"
# Print end status message.
echo
echo "Backup finished"
date
