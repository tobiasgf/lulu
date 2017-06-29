#!/bin/bash
#### Description: concatenate and dereplicate reads from three separate replicates

#### Written by: Tobias Guldberg Frøslev - tobiasgf@bio.ku.dk  on 12-2016

VSEARCH=$(which vsearch)

mkdir single_samples
mv S[0-9][0-9][0-9]* single_samples/

cat single_samples/S031.fas single_samples/S168.fas single_samples/S304.fas  >> S031n.fas
cat single_samples/S043.fas single_samples/S180.fas single_samples/S316.fas  >> S043n.fas
cat single_samples/S128.fas single_samples/S264.fas single_samples/S400.fas  >> S128n.fas
cat single_samples/S093.fas single_samples/S229.fas single_samples/S366.fas  >> S093n.fas
cat single_samples/S047.fas single_samples/S184.fas single_samples/S320.fas  >> S047n.fas
cat single_samples/S048.fas single_samples/S185.fas single_samples/S321.fas  >> S048n.fas
cat single_samples/S049.fas single_samples/S186.fas single_samples/S322.fas  >> S049n.fas
cat single_samples/S050.fas single_samples/S187.fas single_samples/S323.fas  >> S050n.fas
cat single_samples/S051.fas single_samples/S188.fas single_samples/S324.fas  >> S051n.fas
cat single_samples/S052.fas single_samples/S189.fas single_samples/S325.fas  >> S052n.fas
cat single_samples/S053.fas single_samples/S190.fas single_samples/S326.fas  >> S053n.fas
cat single_samples/S055.fas single_samples/S192.fas single_samples/S328.fas  >> S055n.fas
cat single_samples/S056.fas single_samples/S193.fas single_samples/S329.fas  >> S056n.fas
cat single_samples/S057.fas single_samples/S194.fas single_samples/S330.fas  >> S057n.fas
cat single_samples/S058.fas single_samples/S195.fas single_samples/S331.fas  >> S058n.fas
cat single_samples/S090.fas single_samples/S226.fas single_samples/S363.fas  >> S090n.fas
cat single_samples/S059.fas single_samples/S196.fas single_samples/S332.fas  >> S059n.fas
cat single_samples/S060.fas single_samples/S197.fas single_samples/S333.fas  >> S060n.fas
cat single_samples/S061.fas single_samples/S198.fas single_samples/S334.fas  >> S061n.fas
cat single_samples/S062.fas single_samples/S199.fas single_samples/S335.fas  >> S062n.fas
cat single_samples/S063.fas single_samples/S200.fas single_samples/S336.fas  >> S063n.fas
cat single_samples/S064.fas single_samples/S201.fas single_samples/S337.fas  >> S064n.fas
cat single_samples/S065.fas single_samples/S202.fas single_samples/S338.fas  >> S065n.fas
cat single_samples/S066.fas single_samples/S203.fas single_samples/S339.fas  >> S066n.fas
cat single_samples/S072.fas single_samples/S208.fas single_samples/S345.fas  >> S072n.fas
cat single_samples/S073.fas single_samples/S209.fas single_samples/S346.fas  >> S073n.fas
cat single_samples/S074.fas single_samples/S210.fas single_samples/S347.fas  >> S074n.fas
cat single_samples/S075.fas single_samples/S211.fas single_samples/S348.fas  >> S075n.fas
cat single_samples/S076.fas single_samples/S212.fas single_samples/S349.fas  >> S076n.fas
cat single_samples/S077.fas single_samples/S213.fas single_samples/S350.fas  >> S077n.fas
cat single_samples/S110.fas single_samples/S246.fas single_samples/S383.fas  >> S110n.fas
cat single_samples/S111.fas single_samples/S247.fas single_samples/S384.fas  >> S111n.fas
cat single_samples/S112.fas single_samples/S248.fas single_samples/S385.fas  >> S112n.fas
cat single_samples/S113.fas single_samples/S249.fas single_samples/S386.fas  >> S113n.fas
cat single_samples/S114.fas single_samples/S250.fas single_samples/S387.fas  >> S114n.fas
cat single_samples/S115.fas single_samples/S251.fas single_samples/S388.fas  >> S115n.fas
cat single_samples/S116.fas single_samples/S252.fas single_samples/S389.fas  >> S116n.fas
cat single_samples/S117.fas single_samples/S253.fas single_samples/S390.fas  >> S117n.fas
cat single_samples/S118.fas single_samples/S254.fas single_samples/S391.fas  >> S118n.fas
cat single_samples/S119.fas single_samples/S255.fas single_samples/S410.fas  >> S119n.fas
cat single_samples/S120.fas single_samples/S256.fas single_samples/S392.fas  >> S120n.fas
cat single_samples/S121.fas single_samples/S257.fas single_samples/S393.fas  >> S121n.fas
cat single_samples/S122.fas single_samples/S258.fas single_samples/S394.fas  >> S122n.fas
cat single_samples/S123.fas single_samples/S259.fas single_samples/S395.fas  >> S123n.fas
cat single_samples/S124.fas single_samples/S260.fas single_samples/S396.fas  >> S124n.fas
cat single_samples/S134.fas single_samples/S270.fas single_samples/S406.fas  >> S134n.fas
cat single_samples/S125.fas single_samples/S261.fas single_samples/S397.fas  >> S125n.fas
cat single_samples/S126.fas single_samples/S262.fas single_samples/S398.fas  >> S126n.fas
cat single_samples/S127.fas single_samples/S263.fas single_samples/S399.fas  >> S127n.fas
cat single_samples/S129.fas single_samples/S265.fas single_samples/S401.fas  >> S129n.fas
cat single_samples/S130.fas single_samples/S266.fas single_samples/S402.fas  >> S130n.fas
cat single_samples/S131.fas single_samples/S267.fas single_samples/S403.fas  >> S131n.fas
cat single_samples/S132.fas single_samples/S268.fas single_samples/S404.fas  >> S132n.fas
cat single_samples/S135.fas single_samples/S271.fas single_samples/S407.fas  >> S135n.fas
cat single_samples/S136.fas single_samples/S272.fas single_samples/S408.fas  >> S136n.fas
cat single_samples/S137.fas single_samples/S273.fas single_samples/S409.fas  >> S137n.fas
cat single_samples/S023.fas single_samples/S160.fas  >> S023n.fas
cat single_samples/S296.fas  >> S296n.fas
cat single_samples/S054.fas single_samples/S191.fas single_samples/S327.fas  >> S054n.fas
cat single_samples/S105.fas single_samples/S241.fas single_samples/S378.fas  >> S105n.fas
cat single_samples/S016.fas single_samples/S154.fas single_samples/S290.fas  >> S016n.fas
cat single_samples/S017.fas single_samples/S155.fas single_samples/S291.fas  >> S017n.fas
cat single_samples/S018.fas single_samples/S156.fas single_samples/S292.fas  >> S018n.fas
cat single_samples/S069.fas single_samples/S206.fas single_samples/S342.fas  >> S069n.fas
cat single_samples/S070.fas single_samples/S207.fas single_samples/S343.fas  >> S070n.fas
cat single_samples/S019.fas single_samples/S177.fas single_samples/S313.fas  >> S019n.fas
cat single_samples/S020.fas single_samples/S157.fas single_samples/S293.fas  >> S020n.fas
cat single_samples/S021.fas single_samples/S158.fas single_samples/S294.fas  >> S021n.fas
cat single_samples/S022.fas single_samples/S159.fas single_samples/S295.fas  >> S022n.fas
cat single_samples/S024.fas single_samples/S161.fas single_samples/S297.fas  >> S024n.fas
cat single_samples/S009.fas single_samples/S146.fas single_samples/S282.fas  >> S009n.fas
cat single_samples/S010.fas single_samples/S147.fas single_samples/S283.fas  >> S010n.fas
cat single_samples/S011.fas single_samples/S148.fas single_samples/S284.fas  >> S011n.fas
cat single_samples/S012.fas single_samples/S149.fas single_samples/S285.fas  >> S012n.fas
cat single_samples/S013.fas single_samples/S150.fas single_samples/S286.fas  >> S013n.fas
cat single_samples/S014.fas single_samples/S151.fas single_samples/S287.fas  >> S014n.fas
cat single_samples/S040.fas single_samples/S152.fas single_samples/S288.fas  >> S040n.fas
cat single_samples/S068.fas single_samples/S205.fas single_samples/S341.fas  >> S068n.fas
cat single_samples/S015.fas single_samples/S153.fas single_samples/S289.fas  >> S015n.fas
cat single_samples/S001.fas single_samples/S138.fas single_samples/S274.fas  >> S001n.fas
cat single_samples/S002.fas single_samples/S139.fas single_samples/S275.fas  >> S002n.fas
cat single_samples/S003.fas single_samples/S140.fas single_samples/S276.fas  >> S003n.fas
cat single_samples/S004.fas single_samples/S141.fas single_samples/S277.fas  >> S004n.fas
cat single_samples/S005.fas single_samples/S142.fas single_samples/S278.fas  >> S005n.fas
cat single_samples/S006.fas single_samples/S143.fas single_samples/S279.fas  >> S006n.fas
cat single_samples/S007.fas single_samples/S144.fas single_samples/S280.fas  >> S007n.fas
cat single_samples/S008.fas single_samples/S145.fas single_samples/S281.fas  >> S008n.fas
cat single_samples/S067.fas single_samples/S204.fas single_samples/S340.fas  >> S067n.fas
cat single_samples/S101.fas single_samples/S237.fas single_samples/S374.fas  >> S101n.fas
cat single_samples/S102.fas single_samples/S238.fas single_samples/S375.fas  >> S102n.fas
cat single_samples/S103.fas single_samples/S239.fas single_samples/S376.fas  >> S103n.fas
cat single_samples/S104.fas single_samples/S240.fas single_samples/S377.fas  >> S104n.fas
cat single_samples/S106.fas single_samples/S242.fas single_samples/S379.fas  >> S106n.fas
cat single_samples/S107.fas single_samples/S243.fas single_samples/S380.fas  >> S107n.fas
cat single_samples/S108.fas single_samples/S244.fas single_samples/S381.fas  >> S108n.fas
cat single_samples/S109.fas single_samples/S245.fas single_samples/S382.fas  >> S109n.fas
cat single_samples/S133.fas single_samples/S269.fas single_samples/S405.fas  >> S133n.fas
cat single_samples/S078.fas single_samples/S214.fas single_samples/S351.fas  >> S078n.fas
cat single_samples/S091.fas single_samples/S227.fas single_samples/S364.fas  >> S091n.fas
cat single_samples/S079.fas single_samples/S215.fas single_samples/S352.fas  >> S079n.fas
cat single_samples/S080.fas single_samples/S216.fas single_samples/S353.fas  >> S080n.fas
cat single_samples/S081.fas single_samples/S217.fas single_samples/S354.fas  >> S081n.fas
cat single_samples/S082.fas single_samples/S218.fas single_samples/S355.fas  >> S082n.fas
cat single_samples/S083.fas single_samples/S219.fas single_samples/S356.fas  >> S083n.fas
cat single_samples/S084.fas single_samples/S220.fas single_samples/S357.fas  >> S084n.fas
cat single_samples/S085.fas single_samples/S221.fas single_samples/S358.fas  >> S085n.fas
cat single_samples/S092.fas single_samples/S228.fas single_samples/S365.fas  >> S092n.fas
cat single_samples/S094.fas single_samples/S230.fas single_samples/S367.fas  >> S094n.fas
cat single_samples/S095.fas single_samples/S231.fas single_samples/S368.fas  >> S095n.fas
cat single_samples/S096.fas single_samples/S232.fas single_samples/S369.fas  >> S096n.fas
cat single_samples/S097.fas single_samples/S233.fas single_samples/S370.fas  >> S097n.fas
cat single_samples/S098.fas single_samples/S234.fas single_samples/S371.fas  >> S098n.fas
cat single_samples/S099.fas single_samples/S235.fas single_samples/S372.fas  >> S099n.fas
cat single_samples/S100.fas single_samples/S236.fas single_samples/S373.fas  >> S100n.fas
cat single_samples/S039.fas single_samples/S176.fas single_samples/S312.fas  >> S039n.fas
cat single_samples/S086.fas single_samples/S222.fas single_samples/S359.fas  >> S086n.fas
cat single_samples/S087.fas single_samples/S223.fas single_samples/S360.fas  >> S087n.fas
cat single_samples/S088.fas single_samples/S224.fas single_samples/S361.fas  >> S088n.fas
cat single_samples/S089.fas single_samples/S225.fas single_samples/S362.fas  >> S089n.fas
cat single_samples/S044.fas single_samples/S181.fas single_samples/S317.fas  >> S044n.fas
cat single_samples/S071.fas single_samples/S344.fas  >> S071n.fas
cat single_samples/S045.fas single_samples/S182.fas single_samples/S318.fas  >> S045n.fas
cat single_samples/S046.fas single_samples/S183.fas single_samples/S319.fas  >> S046n.fas
cat single_samples/S032.fas single_samples/S169.fas single_samples/S305.fas  >> S032n.fas
cat single_samples/S033.fas single_samples/S170.fas single_samples/S306.fas  >> S033n.fas
cat single_samples/S034.fas single_samples/S171.fas single_samples/S307.fas  >> S034n.fas
cat single_samples/S035.fas single_samples/S172.fas single_samples/S308.fas  >> S035n.fas
cat single_samples/S042.fas single_samples/S179.fas single_samples/S315.fas  >> S042n.fas
cat single_samples/S036.fas single_samples/S173.fas single_samples/S309.fas  >> S036n.fas
cat single_samples/S037.fas single_samples/S174.fas single_samples/S310.fas  >> S037n.fas
cat single_samples/S038.fas single_samples/S175.fas single_samples/S311.fas  >> S038n.fas
cat single_samples/S025.fas single_samples/S162.fas single_samples/S298.fas  >> S025n.fas
cat single_samples/S026.fas single_samples/S163.fas single_samples/S299.fas  >> S026n.fas
cat single_samples/S027.fas single_samples/S164.fas single_samples/S300.fas  >> S027n.fas
cat single_samples/S041.fas single_samples/S178.fas single_samples/S314.fas  >> S041n.fas
cat single_samples/S028.fas single_samples/S165.fas single_samples/S301.fas  >> S028n.fas
cat single_samples/S029.fas single_samples/S166.fas single_samples/S302.fas  >> S029n.fas
cat single_samples/S030.fas single_samples/S167.fas single_samples/S303.fas  >> S030n.fas

for FILE in `ls S[0-9][0-9][0-9]n.fas` ; do
FINAL=${FILE/n.fas/.fas}
   # Dereplicate (vsearch)
    "${VSEARCH}" \
      --derep_fulllength "${FILE}" \
      --sizein \
      --sizeout \
      --fasta_width 0 \
      --relabel_sha1 \
      --output "${FINAL}" > /dev/null
done

rm S*n.fas
