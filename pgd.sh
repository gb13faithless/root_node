##################
# Convert .arp files to .fa
##################
let i = 1

for p in 3seq0/*.arp
do
	java -Xmx1024m -Xms512M -jar PGDSpider2-cli.jar -inputfile "$p" -inputformat ARLEQUIN -outputfile 3seq0/data.fa -outputformat FASTA
	let i++
done