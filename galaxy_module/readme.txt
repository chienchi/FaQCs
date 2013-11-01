This readme file describes the steps you have to perform in order to
install the QC tool within your local installation of Galaxy. 

This operation requires an already running local Galaxy installation (see 
http://g2.bx.psu.edu/ for additional information about galaxy).

Please note that the Galaxy module will provide fewer options than the command
line version because it is meant to be as easy-to-use as possible and, for 
computational reasons.

Below is the step-by-step procedure. Feel free to contact us at 
chienchi@lanl.gov
for any comments or questions.

1) -------

Download and Install R
http://cran.cnr.berkeley.edu/

2) -------

Install the last version of Perl Parallel::ForkManager module 
$ perl -MCPAN -e shell
$ cpan> install Parallel::ForkManager
$ cpan> install String::Approx

3) -------

Copy the Galaxy module xml file (FaQCs.xml) and the QC tool (FaQCs.pl) into 
$ ${galaxy_dir}/tools/FaQCs/
 

4) -------

Add the following three lines 
<section name="FastqQCs" id="FaQCs">
    <tool file="FaQCs/FaQCs.xml"/>
</section>

into
$ ${galaxy_dir}/tool_conf.xml 

for adding the Galaxy left panel entry for FastqQCs section. 

Or you can only add the following line to an existing section.
 <tool file="FaQCs/FaQCs.xml"/>

5) -------

Configure the cluster resource for the tool in a Job runner. 
If your local Galaxy is not hook up with clusters, then ignore this step.

Edit universe_wsgi.ini based on the instruction from 
http://wiki.galaxyproject.org/Admin/Config/Performance/Cluster#Cluster_Resources_Managers

Edit line 44 of the file
$ ${galaxy_dir}/tools/FaQCs/FaQCs.xml

indicating for -t the number of process which should be consistent with the Job runner


6) -------

Restart the Galaxy server



 
