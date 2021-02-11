# ds-mosaic
 Generation of mosaic genome/exome files given two CRAM/BAM genome/exome files.

<b><i>Note:</i></b> 
Dockerfile installs: 
<ul>
<li>bwa version 0.7.17, which is different than publically available bwa internal index files if using bwakit (0.7.12).</li>
<li>samtools 1.11, which may require regeneration of faidx</li>
<li>picard 2.23.8</li>
<li>fastq-tools 0.8.3</li>
</ul>
