##################
# Convert .arp files to .nex
##################

for p in "$dir"/*0
do
	java -Xmx1024m -Xms512M -jar "$PGDdir"/PGDSpider2-cli.jar -inputfile "$p"/*.arp -inputformat ARLEQUIN -outputfile "$p"/data.nex -outputformat NEXUS
done