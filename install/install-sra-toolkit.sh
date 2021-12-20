# Install sra-toolkit.
# This was really annoying to install. The only way to get fasterq-dump to work was by running vdb-config --i at least once
# vdb-config --i causes weird terminal issues. Workaround is to pipe it to null and ignore the exit code =p
set -e

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz
tar -vxzf sratoolkit.2.10.8-ubuntu64.tar.gz
rm sratoolkit.2.10.8-ubuntu64.tar.gz

set +e
echo "q" | ./sratoolkit.2.10.8-ubuntu64/bin/vdb-config -i > /dev/null 2>&1
set -e

ln -s $PWD/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump /usr/local/bin/fasterq-dump
ln -s $PWD/sratoolkit.2.10.8-ubuntu64/bin/prefetch /usr/local/bin/prefetch