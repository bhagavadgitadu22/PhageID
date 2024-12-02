#!/bin/bash

##############################################################################################################
#                                                                                                            #
#  contig_info: a BASH script to estimate standard statistics from FASTA contig files                        #
#                                                                                                            #
   COPYRIGHT="Copyright (C) 2018-2021 Institut Pasteur"                                                      #
#                                                                                                            #
#  This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU  #
#  General Public License as published by the Free Software Foundation, either version 3 of the License, or  #
#  (at your option) any later version.                                                                       #
#                                                                                                            #
#  This program is distributed in the hope that it will be useful,  but WITHOUT ANY WARRANTY;  without even  #
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public  #
#  License for more details.                                                                                 #
#                                                                                                            #
#  You should have received a copy of the  GNU General Public License along with this program.  If not, see  #
#  <http://www.gnu.org/licenses/>.                                                                           #
#                                                                                                            #
#  Contact:                                                                                                  #
#   Alexis Criscuolo                                                            alexis.criscuolo@pasteur.fr  #
#   Genome Informatics & Phylogenetics (GIPhy)                                             giphy.pasteur.fr  #
#   Bioinformatics and Biostatistics Hub      research.pasteur.fr/team/bioinformatics-and-biostatistics-hub  #
#   Dpt. Biologie Computationnelle                     research.pasteur.fr/department/computational-biology  #
#   Institut Pasteur, Paris, FRANCE                                                     research.pasteur.fr  #
#                                                                                                            #
##############################################################################################################

##############################################################################################################
#                                                                                                            #
# ============                                                                                               #
# = VERSIONS =                                                                                               #
# ============                                                                                               #
#                                                                                                            #
  VERSION=2.01.210324ac                                                                                      #
# + modified temp file creation using TMPDIR variable                                                        #
# + adding option -h                                                                                         #
#                                                                                                            #
# VERSION=2.0.210315ac                                                                                       #
# + adding option -r                                                                                         #
# + adding trap when interrupting script                                                                     #
#                                                                                                            #
# VERSION=1.1.201007ac                                                                                       #
# + estimating auN (also called E-size)                                                                      #
#                                                                                                            #
# VERSION=1.0.190426ac                                                                                       #
# + options -l and -d (i.e. printing sequence lengths and length distribution) are no longer supported       #
# + residue count always computed (option -r discarded)                                                      #
# + ultrafast residue count (based on tr + wc)                                                               #
# + estimating %AT, %GC, L50, L75, L90                                                                       #
# + faster estimation of the sequence length statistics (100% awk)                                           #
# + ability to read multiple input files                                                                     #
#                                                                                                            #
# VERSION=0.3.180515ac                                                                                       #
#                                                                                                            #
##############################################################################################################
  
##############################################################################################################
#                                                                                                            #
#  ================                                                                                          #
#  = INSTALLATION =                                                                                          #
#  ================                                                                                          #
#                                                                                                            #
#  Just give the execute permission to the script contig_info.sh with the following command line:            #
#                                                                                                            #
#   chmod +x contig_info.sh                                                                                  #
#                                                                                                            #
##############################################################################################################

#############################################################################################################
#                                                                                                           #
#  ================                                                                                         #
#  = FUNCTIONS    =                                                                                         #
#  ================                                                                                         #
#                                                                                                           #
# = mandoc() =============================================================================================  #
#   prints the doc                                                                                          #
#                                                                                                           #
mandoc() {
  echo -e "\n\033[1m contig_info v$VERSION           $COPYRIGHT\033[0m";
  cat <<EOF

 USAGE:  contig_info.sh  [options]  <contig_files> 

  where 'options' are:

   -m <int>    minimum contig length; every contig sequence of length shorter
               than this cutoff will be discarded (default: 1)
   -g <int>    expected genome size  for computing auNG  and {N,L}G{50,75,90}
               values instead of auN and {N,L}{50,75,90} ones, respectively
   -r          residue content statistics for each contig sequence instead of 
               global statistics
   -t          tab-delimited output
   -h          prints this help and exits

EOF
}
#                                                                                                           #
# = randomfile() =========================================================================================  #
#   creates a random file within $TEMP_DIR and returns its name                                             #
#                                                                                                           #
randomfile() {
  echo $(mktemp -t -p ${TMPDIR:-/tmp}) ;
}
#                                                                                                           #
#############################################################################################################

#############################################################################################################
####                                                                                                     ####
#### INITIALIZING PARAMETERS AND READING OPTIONS                                                         ####
####                                                                                                     ####
#############################################################################################################
if [ $# -lt 1 ]; then mandoc ; exit 1 ; fi

export LC_ALL=C;

MIN_CONTIG_LGT=1;
GENOME_SIZE=0;
RES_CONTENT=false;
TSVOUT=false;

while getopts m:g:trh option
do
  case $option in
    m)  MIN_CONTIG_LGT=$OPTARG ;;
    g)  GENOME_SIZE=$OPTARG    ;;
    r)  RES_CONTENT=true       ;;
    t)  TSVOUT=true            ;;
    h)  mandoc ;        exit 0 ;;
    \?) mandoc ;        exit 1 ;;
  esac
done
if ! [[ $MIN_CONTIG_LGT =~ ^[0-9]+$ ]]; then echo " incorrect value (option -m): $MIN_CONTIG_LGT"                           >&2 ; exit 1 ; fi
if [ $MIN_CONTIG_LGT -lt 1 ];           then echo " the min contig length threshold must be a positive integer (option -m)" >&2 ; exit 1 ; fi
if ! [[ $GENOME_SIZE =~ ^[0-9]+$ ]];    then echo " incorrect value (option -g): $GENOME_SIZE"                              >&2 ; exit 1 ; fi
if [ $GENOME_SIZE -lt 0 ];              then echo " the expected genome size must be a positive integer (option -g)"        >&2 ; exit 1 ; fi


#############################################################################################################
####                                                                                                     ####
#### CONTIG INFO                                                                                         ####
####                                                                                                     ####
#############################################################################################################
if $TSVOUT
then
  if ! $RES_CONTENT
  then
    CSVCAPT="#File\tNseq\tNres\tA\tC\tG\tT\tN\t%A\t%C\t%G\t%T\t%N\t%AT\t%GC\tMin\tQ25\tMed\tQ75\tMax\tAvg\tauN\tN50\tN75\tN90\tL50\tL75\tL90";
    [ $GENOME_SIZE -ne 0 ] && CSVCAPT="$CSVCAPT\tExpSize";
  else
    CSVCAPT="#File\tSeq\tNres\tA\tC\tG\tT\tN\t%A\t%C\t%G\t%T\t%N\t%AT\t%GC\tPval";
  fi
  echo -e "$CSVCAPT" ;
fi

SEQS=$(randomfile); ## txt file with one sequence per line
trap "echo -n interrupting ... ; rm -f $SEQS &>/dev/null ; echo ; exit " SIGINT

for INFILE in "$@"
do
  if [ ! -e $INFILE ] || [ ! -f $INFILE ] || [ ! -r $INFILE ]; then continue; fi
    
  tr -d '\15\32' < $INFILE | awk -v thr=$MIN_CONTIG_LGT '/^>/{if(length(s)>=thr)print s;s="";next}{s=s$0}END{if(length(s)>=thr)print s}' | tr '[:lower:]' '[:upper:]' > $SEQS ;
  
  A=$(tr -cd A < $SEQS | wc -c);
  C=$(tr -cd C < $SEQS | wc -c); 
  G=$(tr -cd G < $SEQS | wc -c);
  T=$(tr -cd T < $SEQS | wc -c); 
  ACGT=$(( $A + $C + $G + $T ));

  if ! $RES_CONTENT
  then
    R=$(tr -d [:cntrl:] < $SEQS | wc -c);
    S=$(wc -l < $SEQS);
    AVG=$(bc -l<<<"scale=2;$R/$S" | sed 's/^\./0./');
    fA=$(bc -l <<<"scale=2;100*$A/$R" | sed 's/^0$/0.00/;s/^\./0./');
    fC=$(bc -l <<<"scale=2;100*$C/$R" | sed 's/^0$/0.00/;s/^\./0./');
    fG=$(bc -l <<<"scale=2;100*$G/$R" | sed 's/^0$/0.00/;s/^\./0./');
    fT=$(bc -l <<<"scale=2;100*$T/$R" | sed 's/^0$/0.00/;s/^\./0./');
    N=0; [ $R -ne $ACGT ] && N=$(tr -cd N < $SEQS | wc -c);
    fN=$(bc -l <<<"scale=2;100*$N/$R" | sed 's/^0$/0.00/;s/^\./0./');
    fGC=$(bc -l <<<"scale=2;100*($C+$G)/$ACGT" | sed 's/^0$/0.00/;s/^\./0./');
    fAT=$(bc -l <<<"scale=2;100-$fGC" | sed 's/^0$/0.00/;s/^\./0./');
    ER=$R; [ $GENOME_SIZE != 0 ] && ER=$GENOME_SIZE;
    STATS=$(awk '{print length}' $SEQS | sort -rn | awk -v g=$ER '{l[++n]=$0;aun+=$0*$0}
                                                                  END{OFMT="%f";g50=g/2;g75=3*g/4;g90=9*g/10;i=s=n50=n75=n90=0;
                                                                      while(++i<=n&&n90==0){s+=l[i];n50==0&&s>=g50&&n50=l[i]+(l50=i);n75==0&&s>=g75&&n75=l[i]+(l75=i);n90==0&&s>=g90&&n90=l[i]+(l90=i)}
                                                                      n90==0&&n90=l[n]+(l90=n);n75==0&&n75=l[n]+(l75=n);n50==0&&n50=l[n]+(l50=n);
                                                                      iq3=int(n/4+1);iq2=int(n/2+1);iq1=int(3*n/4+1);
                                                                      print(n50-l50)"\t"(n75-l75)"\t"(n90-l90)"\t"l50"\t"l75"\t"l90"\t"l[1]"\t"l[iq3]"\t"l[iq2]"\t"l[iq1]"\t"l[n]"\t"int(0.5+aun/g)}');
    N50=$(cut -f1  <<<"$STATS"); N75=$(cut -f2 <<<"$STATS");  N90=$(cut -f3  <<<"$STATS"); 
    L50=$(cut -f4  <<<"$STATS"); L75=$(cut -f5 <<<"$STATS");  L90=$(cut -f6  <<<"$STATS"); 
    Q75=$(cut -f8  <<<"$STATS"); Q50=$(cut -f9 <<<"$STATS");  Q25=$(cut -f10 <<<"$STATS"); 
    MIN=$(cut -f11 <<<"$STATS"); MAX=$(cut -f7 <<<"$STATS");  AUN=$(cut -f12 <<<"$STATS");

    if ! $TSVOUT
    then
      echo ;
      echo "File                           $(basename $INFILE)" ;
      echo ;
      echo "Number of sequences            $S" ;
      echo ;
      echo "Residue counts:" ;
      echo "  Number of A's                $A  $fA %" ;
      echo "  Number of C's                $C  $fC %" ;
      echo "  Number of G's                $G  $fG %" ;
      echo "  Number of T's                $T  $fT %" ;
      echo "  Number of N's                $N  $fN %" ;
      echo "  Total                        $R" ;
      echo ;
      echo "  %AT                          $fAT %" ;
      echo "  %GC                          $fGC %" ;
      echo ;
      echo "Sequence lengths:" ;
      echo "  Minimum                      $MIN" ;
      echo "  Quartile 25%                 $Q25" ;
      echo "  Median                       $Q50" ;
      echo "  Quartile 75%                 $Q75" ;
      echo "  Maximum                      $MAX" ;
      echo "  Average                      $AVG" ;
      echo ;
      echo "Contiguity statistics:" ;
      echo "  auN                          $AUN" ;
      echo "  N50                          $N50" ;
      echo "  N75                          $N75" ;
      echo "  N90                          $N90" ;
      echo "  L50                          $L50" ;
      echo "  L75                          $L75" ;
      echo "  L90                          $L90" ;
      if [ $GENOME_SIZE -ne 0 ]; then echo "  Expected genome size         $GENOME_SIZE"; fi
      echo ;
    else
      CSVLINE="$(basename $INFILE)\t$S\t$R\t$A\t$C\t$G\t$T\t$N\t$fA%\t$fC%\t$fG%\t$fT%\t$fN%\t$fAT%\t$fGC%\t$MIN\t$Q25\t$Q50\t$Q75\t$MAX\t$AVG\t$AUN\t$N50\t$N75\t$N90\t$L50\t$L75\t$L90";
      [ $GENOME_SIZE -ne 0 ]&&CSVLINE="$CSVLINE\t$ER";
      echo -e "$CSVLINE" ;
    fi

    continue ;
  fi

  
  ## residue details

  WS=200;              # window size
  NS=5000;             # no. samples from all the contigs
  GCALL=$(randomfile); # %GC estimated from $NS segments of length $WS regularly sampled from all the contigs
  tr ACGT WSSW < $SEQS | fold -w $WS | tr -cd WS'\n' | awk -v ws=$WS -v step=$(( 1 + $ACGT / ($WS*$NS) )) '(length()!=ws){next}(++n==step){printf("%.6f\n",gsub("S","S")/ws);n=0}' > $GCALL ;
  NS=500;              # no. samples from each contig
  GCSEQ=$(randomfile); # %GC estimated from $NS segments of length $WS regularly sampled from each contig
  
  trap "echo -n interrupting ... ; rm -f $SEQS $GCALL $GCSEQ &>/dev/null ; echo ; exit " SIGINT

  tr -d '\15\32' < $INFILE | awk '!/^>/{s=s$0;next}(s!=""){print toupper(s);s=""}{print}END{print toupper(s)}' | paste - - > $SEQS ; ## tsv file: FASTA header \t contig

  while IFS=$'\t' read -r fh seq
  do
    R=${#seq};
    if [ $R -lt $MIN_CONTIG_LGT ]; then continue; fi

    NAME=$(tr -d '>' <<<"$fh");
    A=$(tr -cd A <<<"$seq" | wc -c); fA=$(bc -l <<<"scale=2;100*$A/$R" | sed 's/^0$/0.00/;s/^\./0./');
    C=$(tr -cd C <<<"$seq" | wc -c); fC=$(bc -l <<<"scale=2;100*$C/$R" | sed 's/^0$/0.00/;s/^\./0./');
    G=$(tr -cd G <<<"$seq" | wc -c); fG=$(bc -l <<<"scale=2;100*$G/$R" | sed 's/^0$/0.00/;s/^\./0./');
    T=$(tr -cd T <<<"$seq" | wc -c); fT=$(bc -l <<<"scale=2;100*$T/$R" | sed 's/^0$/0.00/;s/^\./0./');
    ACGT=$(( $A + $C + $G + $T ));
    N=0; [ $R -ne $ACGT ] && N=$(tr -cd N <<<"$seq" | wc -c);
    fN=$(bc -l <<<"scale=2;100*$N/$R" | sed 's/^0$/0.00/;s/^\./0./');
    fGC=$(bc -l <<<"scale=2;100*($C+$G)/$ACGT" | sed 's/^0$/0.00/;s/^\./0./');
    fAT=$(bc -l <<<"scale=2;100-$fGC" | sed 's/^0$/0.00/;s/^\./0./');

    ### 1. getting (up to) $NS segments of length $WS regularly sampled from $seq and saving all the corresponding %GC in $GCSEQ
    tr ACGT WSSW <<<"$seq" | fold -w $WS | tr -cd WS'\n' | awk -v ws=$WS -v step=$(( 1 + ($R / ($WS*$NS)) )) '(length()!=ws){next}(++n==step){printf("%.6f\n",gsub("S","S")/ws);n=0}' > $GCSEQ ; 
    ### 2. comparing the %GC values sampled from $seq (file $GCSEQ) to the %GC values sampled from all the contigs (file $GCALL)
    ###    the adequation between the two sets of %GC values is assessed by a Mann-Whitney U test
    pv=$(awk 'function abs(x){return (x<0)?-x:x}
              function erf(x){q=t=1.0/(1+0.47047*x);sum=0.3480242*t;sum-=0.0958798*(t*=q);sum+=0.7478556*(t*=q);return 1.0-(sum*exp(-x*x))}
              (FNR==NR){a1[++n1]=$0;next}{a2[++n2]=$0}
              END{while(++i2<=n2){v2=a2[i2];i1=0;while(++i1<=n1)(v2>a1[i1])&&++u}
                  if(n2==0){print"n/a";exit}
                  if(n2==1){p=u/n1;p=abs(0.5-p);printf("%.4f",1-2*p);exit}
                  mean=n1*n2/2;(u>mean)&&u=2*mean-u;sd=sqrt(mean*(n1+n2+1)/6);z=abs(u-mean)/sd;p=(1+erf(z/sqrt(2)))/2;p=2*(1-p); printf("%.4f",p)}' $GCALL $GCSEQ);
    false && awk 'function abs(x){return (x<0)?-x:x}
                  function erf(x){q=t=1.0/(1+0.47047*x);sum=0.3480242*t;sum-=0.0958798*(t*=q);sum+=0.7478556*(t*=q);return 1.0-(sum*exp(-x*x))}
                  (FNR==NR){a1[++n1]=$0;next}{a2[++n2]=$0}
                  END{while(++i2<=n2){v2=a2[i2];i1=0;while(++i1<=n1)(v2>a1[i1])&&++u}
                      if(n2==0){print"n/a";exit}
                      if(n2==1){p=u/n1;p=abs(0.5-p);print n1" "n2" "u" "(u/n1)" "(1-2*p);exit}
                      mean=n1*n2/2;(u>mean)&&u=2*mean-u;sd=sqrt(mean*(n1+n2+1)/6);z=abs(u-mean)/sd;p=(1+erf(z/sqrt(2)))/2;p=2*(1-p); print n1" "n2" "u" "mean" "sd" "z" "p}' $GCALL $GCSEQ
    
    if ! $TSVOUT
    then
      echo ;
      echo "File                           $(basename $INFILE)" ;
      echo ;
      echo "Sequence                       $NAME" ;
      echo ;
      echo "Residue counts:" ;
      echo "  Number of A's                $A  $fA %" ;
      echo "  Number of C's                $C  $fC %" ;
      echo "  Number of G's                $G  $fG %" ;
      echo "  Number of T's                $T  $fT %" ;
      echo "  Number of N's                $N  $fN %" ;
      echo "  Total                        $R" ;
      echo ;
      echo "  %AT                          $fAT %" ;
      echo "  %GC                          $fGC %" ;
      echo ;
      echo "Composition test p-value:      $pv" ;
      echo ;
    else
      CSVLINE="$(basename $INFILE)\t$NAME\t$R\t$A\t$C\t$G\t$T\t$N\t$fA%\t$fC%\t$fG%\t$fT%\t$fN%\t$fAT%\t$fGC%\t$pv";
      echo -e "$CSVLINE" ;
    fi
  done < $SEQS
  
done

rm -f $SEQS $GCALL $GCSEQ ;

exit ;

  
 
