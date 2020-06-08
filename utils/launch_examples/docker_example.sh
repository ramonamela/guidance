#!/bin/bash

export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/
export MPI_HOME=/usr/lib/openmpi
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/include/:${LD_LIBRARY_PATH}
cat /etc/environment
source /etc/profile.d/compss.sh
cat /etc/profile.d/compss.sh

amount_cpus=$(cat /proc/cpuinfo | grep processor | wc | awk '{ print $1 }')

sed -i "s@<ComputingUnits>8</ComputingUnits>@<ComputingUnits>${amount_cpus}</ComputingUnits>@g" /guidance/utils/xml_files/resources.xml


source /guidance/utils/set_environment.sh

#apt-get install -y net-tools

mkdir -p ~/.ssh
ssh-keygen -t rsa -f "~/.ssh/id_rsa"
cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
chmod 700 ~/.ssh
chmod 600 ~/.ssh/authorized_keys
sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config
sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
sudo /usr/sbin/sshd

#echo "export SHAPEITBINARY=\"/usr/bin/shapeit\"" >> ~/.bashrc
#echo "export PLINKBINARY=\"/usr/bin/plink\"" >> ~/.bashrc
#echo "export QCTOOLBINARY=\"/usr/bin/qctool1.4\"" >> ~/.bashrc
#echo "export EAGLEBINARY=\"/usr/bin/eagle\"" >> ~/.bashrc
#echo "export IMPUTE2BINARY=\"/usr/bin/impute2\"" >> ~/.bashrc
#echo "export QCTOOLBINARY=\"/usr/bin/qctool1.4\"" >> ~/.bashrc
#echo "export SNPTESTBINARY=\"/usr/bin/snptest_v2.5\"" >> ~/.bashrc
#echo "export MINIMAC3BINARY=\"/usr/bin/minimac3\"" >> ~/.bashrc
#echo "export MINIMAC4BINARY=\"/usr/bin/minimac4\"" >> ~/.bashrc
#echo "export TABIXBINARY=\"/usr/bin/tabix\"" >> ~/.bashrc
#echo "export BGZIPBINARY=\"/usr/bin/bgzip\"" >> ~/.bashrc
#echo "export BCFTOOLSBINARY=\"/usr/bin/bcftools\"" >> ~/.bashrc
#echo "export SAMTOOLSBINARY=\"/usr/bin/samtools\"" >> ~/.bashrc
#echo "export LC_ALL=\"C\"" >> ~/.bashrc
#echo "export RSCRIPTDIR=\"/TOOLS/R_scripts/\"" >> ~/.bashrc
#echo "export RSCRIPTBINDIR=\"/usr/bin/\"" >> ~/.bashrc
#
#echo "export SHAPEITBINARY=\"/usr/bin/shapeit\"" >> /etc/environment
#echo "export PLINKBINARY=\"/usr/bin/plink\"" >> /etc/environment
#echo "export QCTOOLBINARY=\"/usr/bin/qctool1.4\"" >> /etc/environment
#echo "export EAGLEBINARY=\"/usr/bin/eagle\"" >> /etc/environment
#echo "export IMPUTE2BINARY=\"/usr/bin/impute2\"" >> /etc/environment
#echo "export QCTOOLBINARY=\"/usr/bin/qctool1.4\"" >> /etc/environment
#echo "export SNPTESTBINARY=\"/usr/bin/snptest_v2.5\"" >> /etc/environment
#echo "export MINIMAC3BINARY=\"/usr/bin/minimac3\"" >> /etc/environment
#echo "export MINIMAC4BINARY=\"/usr/bin/minimac4\"" >> /etc/environment
#echo "export TABIXBINARY=\"/usr/bin/tabix\"" >> /etc/environment
#echo "export BGZIPBINARY=\"/usr/bin/bgzip\"" >> /etc/environment
#echo "export BCFTOOLSBINARY=\"/usr/bin/bcftools\"" >> /etc/environment
#echo "export SAMTOOLSBINARY=\"/usr/bin/samtools\"" >> /etc/environment
#echo "export LC_ALL=\"C\"" >> /etc/environment
#echo "export RSCRIPTDIR=\"/TOOLS/R_scripts/\"" >> /etc/environment
#echo "export RSCRIPTBINDIR=\"/usr/bin/\"" >> /etc/environment

export SHAPEITBINARY="/usr/bin/shapeit"
export PLINKBINARY="/usr/bin/plink"
export QCTOOLBINARY="/usr/bin/qctool1.4"
export EAGLEBINARY="/usr/bin/eagle"
export IMPUTE2BINARY="/usr/bin/impute2"
export QCTOOLBINARY="/usr/bin/qctool1.4"
export SNPTESTBINARY="/usr/bin/snptest_v2.5"
export MINIMAC3BINARY="/usr/bin/minimac3"
export MINIMAC4BINARY="/usr/bin/minimac4"
export TABIXBINARY="/usr/bin/tabix"
export BGZIPBINARY="/usr/bin/bgzip"
export BCFTOOLSBINARY="/usr/bin/bcftools"
export SAMTOOLSBINARY="/usr/bin/samtools"
export LC_ALL="C"
export RSCRIPTDIR="/TOOLS/R_scripts/"
export RSCRIPTBINDIR="/usr/bin/"

pushd /guidance
mvn install
popd

runcompss \
	--project="/guidance/utils/xml_files/project.xml" \
  --resources="/guidance/utils/xml_files/resources.xml" \
  --classpath="/guidance/guidance.jar" \
	--jvm_workers_opts="-Dcompss.worker.removeWD=false" \
  guidance.Guidance -config_file /guidance/utils/conf_examples/config_base_example.file


echo "Done!"

sleep 50000
