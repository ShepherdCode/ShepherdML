#!/bin/bash
#read.file.line.by.line.sh

while read line
do
#echo $line | md5sum | awk '{print $1}'
echo $line | md5
done

exit

# This shows we have 2 sequences that occur 7 times
$ grep -v '^>' pcRNA.fasta | ./hash_by_line.sh | sort | uniq -c | sort -n | tail
   2 fa3f53471c7a1c1e2412b563ae4e35b8
   2 fefa3ddfaacd27e593ce6f56db0e07f5
   3 1e89cdd414f815a3eb05bc1f3efaf65c
   3 28ec4da6b0e7224769b5d56c97b0184e
   3 4e3b0a89de16701eae81b517af3d62c8
   3 56ce36e6ce17df097fee817f25b5212c
   3 65f106060dac7f0e39dc10144f07137f
   5 b0318e1431c176121fbb19e37d3f5164
   7 2b0e28c878799780914f11d90f2294e9
   7 ffd824f50452aa9e3421b6c0d049b626
