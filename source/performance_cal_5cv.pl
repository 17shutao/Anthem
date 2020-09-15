use 5.010;
use IO::File;
use Cwd;
#$arg1=shift;
#$arg2=shift;
#print "$arg1, $arg2";
	&performance_cal_by_perl_combine();

	sub performance_cal_by_perl_combine(){
		my $path=getcwd();
		&performance_cal_by_cutoff($path); #step 1  #通过变换cutoff,来生成用于perl计算auc的文件
		my $min=&select_cutoff($path);  	# step 2
		&svm_mcc_calculate_by_regression($min,$path); #step 3
	}
	
	
	
	sub select_cutoff(){
		my $path=shift;
		my $in=new IO::File("<$path/myresult_AUC.txt") or die; 
  		my $out=new IO::File(">$path/min.txt") or die; 
  		my @lines=<$in>;
  		my $min=2;
  		my %hs_d;
  		for my $l (@lines){
  			my @data=split/\s+/,$l;
  			my $d=sqrt(($data[1]-0)**2+($data[2]-1)**2);
  			$hs_d{$data[0]}=$d;
  		}
  		for my $k (sort{$hs_d{$a} <=> $hs_d{$b} } (keys %hs_d)){
			print $out "$k $hs_d{$k}\n";  
			return $k;			
  		}
	}
		
	sub performance_cal_by_cutoff(){  #通过变换cutoff,来生成用于perl计算auc的文件
        my $path=shift;
  			my $in1=new IO::File("<$path/result_regression.txt") or die;   # change here
  			my $out1=new IO::File(">$path/myresult_AUC.txt") or die;  # change here
  			my $out2=new IO::File(">$path/myresult_RP.txt") or die; # change here
  			for(my $y=-1;$y<=1;$y+=0.001){
          my $tp=0;
    			my $fp=0;
    			my $tn=0;
    			my $fn=0;
    			my $fpr=0;
    			my $recall=0;
    			my $precision=0;
         	while(defined(my $line1=<$in1>)) {  
         		chomp $line1;
         		my @a=split/\s+/,$line1;
         		#print "$a[0] $a[1]\n";
         		if(($a[0] == 1)&&($a[1] >= $y)){
          		$tp++;
          	}
          	if(($a[0] == 1)&&($a[1] < $y)){
          		$fn++;
          	}
          	if(($a[0] == -1)&&($a[1] >= $y)){
          		$fp++;
          	}
          	if(($a[0] == -1)&&($a[1] < $y)){
          		$tn++;
          	}
         	}
         	seek $in1,0,0;
         	print "$tp $fn $fp $tn\n";
    			my $sum=$tp+$tn+$fp+$fn;
    			my $acc=($tp+$tn)/($tp+$tn+$fp+$fn);    
          my $sensitivity=$tp/($tp+$fn);
          my $specificity=$tn/($tn+$fp);
          my $mcc;
          if(($tp eq 0)&&($fp eq 0)){
            $precision=1;
            #$mcc=0;
          }
          else{
            $precision=$tp/($tp+$fp);
            #$mcc=($tp*$tn-$fp*$fn)/sqrt(($tn+$fn)*($tn+$fp)*($tp+$fn)*($tp+$fp));
          }
          $recall=$tp/($tp+$fn);#Tpr
          $fpr=$fp/($fp+$tn);
          print $out1 "$y $fpr $recall  $tp  $tn $fp  $fn    $sum\n";
          print $out2 "$y $recall  $precision $tp  $tn $fp  $fn   $sum  $specificity\n";
		  }
	}
	
	
	
	sub svm_mcc_calculate_by_regression(){
		my $cutoff=shift||0;
		my $path=shift;
    		my $out=new IO::File(">$path/performance.txt") or die "$!";
    		my ($result,$precision,$acc,$sensitivity,$specificity,$total,$mcc,$tp,$fp,$tn,$fn,);
            		$tp=0;
            		$fp=0;
            		$tn=0;
            		$fn=0;	
 				my $in_test=new IO::File("<$path/result_regression.txt") or die "$!";#change the name
                		while(defined(my $l_test=<$in_test>)){
        	    			my @list_test=split/\s+/,$l_test;
        	    			my $gold=$list_test[0];
                			my $result=$list_test[1];
					if(($result >= $cutoff)&&($gold == 1)){
						$tp++;
					}
					elsif(($result < $cutoff)&&($gold == 1)){
						$fn++;
					}
                			elsif(($result < $cutoff)&&($gold == -1)){
        						$tn++;
        					}
					elsif(($result >= $cutoff)&&($gold == -1)){
						$fp++;
					}
					else{
						print "error!";
						die;
					}
				}

			$mcc=($tp*$tn-$fp*$fn)/sqrt(($tn+$fn)*($tn+$fp)*($tp+$fn)*($tp+$fp));
			$acc=($tp+$tn)/($tp+$tn+$fp+$fn);    
			$sensitivity=$tp/($tp+$fn);
			$specificity=$tn/($tn+$fp);
			$precision=$tp/($tp+$fp);
       			# print $out "$k times mcc_mean is $mcc_mean\n";
       	 		print $out "mcc acc sensitivity specificity precision tp fp tn fn\n";
        		print $out "$mcc $acc $sensitivity $specificity $precision $tp $fp $tn $fn\n";
       }
	
	