# Nanopore long-read の端に見える２塩基の繰り返し（poly-A/T/G/C含む）を消す。
my $pol_len = 10;		# 最小のリピート回数。小さいと正常な配列も取っちゃう。８とした場合はAGAGAGAGAGAGAGAG (16b) より端を削る
my $shortlimit = 100;	# 短すぎるreadの判定基準。これより短くなったやつは消す。長いやつを処理する。
my $target_unit = 150;	# 端から何塩基単位に判定するか。長いと雑。短いと丁寧だが時間がかかる
my $trimout_file = 'endtrim_trimawayseq.txt';
my $report_file = 'endtrim_report.txt';

my $pass = $ARGV[0];
my $outpass = $ARGV[1];
if($pass eq ''){
	print 'Usage'."\n";
	print 'perl endtrim.pl [input_directry] [output_directry]'."\n";
	print 'if [output directry] is not provided, ./endtrim_out/ directry may be created'."\n";
	print 'only fastq data is applicable (not .gz)'."\n";
	exit(0);
}
if($outpass eq ''){$outpass = './endtrim_out/';}

print $pass."\n";
print $outpass."\n";

$|=1;

my @nuc;			# 塩基の記号。４種類
$nuc[0] = 'A';
$nuc[1] = 'T';
$nuc[2] = 'G';
$nuc[3] = 'C';
my $nucs=4;

my @polyn;			# 判定用の文字列の配列を生成する。指定回数の２塩基繰り返し。
for(my $k=0;$k<$nucs;$k++){
for(my $j=0;$j<$nucs;$j++){
	$polyn[$k*$nucs+$j] = '';
for(my $i=0;$i<$pol_len;$i++){
	$polyn[$k*$nucs+$j] = $polyn[$k*$nucs+$j].$nuc[$k].$nuc[$j];	# とりあえず同じ塩基の繰り返し
}}}
my $kinds = $nucs * $nucs;	# ２塩基の繰り返しの種類。４種類x４種類=16種類

my @file;			# 対象ファイルを探す。
my @outfile;
my $files=0;
opendir DIR, $pass;
my @dir = readdir(DIR);
foreach my $f(@dir){
	if($f =~ /endtrim/){}else{					# コレで出力されるやつは外す
	if($f =~ /\.fastq$/){
		$f =~ s/\.fastq$//;						# .fastq で終わるファイルを対象とする
		$file[$files]=$pass.$f.'.fastq';
		$outfile[$files]=$outpass.$f.'_endtrim.fastq';
		$files++;
	}}
}
close DIR;

# 確認
for(my $f=0;$f<$files;$f++){
	print $file[$f]."\n";
	print $outfile[$f]."\n";
}
print $trimout_file."\n";
print $report_file."\n";
#exit(0);

open TRM, '>'.$trimout_file;

my $seq;	# 引数渡しが面倒くさいのでグローバル
my $qual;
sub edge(){								# トリミングする機能は関数化しておく
	if(length($seq) < $shortlimit){return(0);}
	my $trim = 0;	# トリミングしたかどうかのフラグ。returnで返す。
	
	my $flag = 1;				# 左側をトリミング
	while($flag==1){	# 長さに合わせて柔軟にトリミングするために短い判定を繰り返す
		my $pos=(-1);							# トリミング位置。初期値は(-1)。rindexでマッチなしのときの返り値
		my $str = substr($seq,0,$target_unit);	# 判定する文字列。左端から指定の塩基長分
		for(my $j=0;$j<$kinds;$j++){			# 判定用配列と一致するもののうち一番右のもの（位置の値が最大値）を探す
			my $posi = rindex($str, $polyn[$j]);	# rindex関数を使う。ヒットしたうち一番右の位置を返す。ヒットしないときは(-1)を返す
			if($pos < $posi){$pos = $posi;}
		}
		if($pos >= 0){			# トリミング位置がみつかったときは >=0 の値を取る
			$pos = $pos + $pol_len*2;	# トリミング位置はヒット位置の値よりも検索サイズ分右にずれた位置
			if($pos > $target_unit - $pol_len*2){$pos=$pos-$pol_len*2;}		# ウインドウの境目に当たったとき、次のラウンドで残さずトリムするための処理
			print TRM 'left: ('.$pos.' b)'.substr($seq, 0, $pos)."\n";
			$seq = substr($seq, $pos);		# 第３引数を省略すると最後まで。
			$qual= substr($qual,$pos);
			$trim=1;
			$trimbases += $pos;
			if(length($seq) < $shortlimit){return(0);}		# 削って小さくなったら無かったことにするのでゼロを返す
		}else{$flag=0;}
	}
	if($trim==1){print TRM 'left: '.substr($seq, 0, 100).'...'."\n\n";}
	
	my $trim_=0;
	$flag = 1;					# 右側をトリミング
	while($flag==1){
		$pos = 0;			# トリミング位置の初期値。index()がヒットなし場合(-1)を返すが、これが最もトリム帳が長くなる一番左を意味してしまうのでややこしい
		$str = substr($seq,length($seq)-$target_unit);		# 判定する文字列。右側。
		for(my $j=0;$j<$kinds;$j++){
			my $posi = index($str, $polyn[$j]);		# index関数を使う。ヒットしたうち一番左の位置を返す。ヒットしないとさらに１個左(-1)を返す。
			if($posi > (-1)){						# index関数の返り値を判定する必要がある。0は最も多くトリムする位置で、(-1)はトリムしない判定。
				$posi = $target_unit - $posi;	# 反転する。トリムする塩基長を意味する値に変わる。
				if($pos < $posi){$pos=$posi;}		# トリムする塩基長が最大の値を探す。
			}
		}
		if($pos > 0){	# ヒットしていれば反転されているのでトリム長を意味する。判定用配列の長さより短いトリム長にはならない。ヒットしなかった場合は０のまま。
			if($pos > $target_unit - $pol_len*2){$pos=$pos-$pol_len*2;}		# 境目で残っちゃう対策
			print TRM 'right:('.$pos.' b)'.substr($seq, length($seq)-$pos)."\n";
			$seq = substr($seq, 0, length($seq )-$pos);
			$qual= substr($qual,0, length($qual)-$pos);		# 左端からトリム長分減った位置までにする
			$trim=1;
			$trim_=1;
			$trimbases += $pos;
			if(length($seq) < $shortlimit){return(0);}		# 削って小さくなったら無かったことにするのでゼロを返す
		}else{$flag=0;}
	}
	if($trim_==1){print TRM 'right: ....'.substr($seq, length($seq)-100)."\n\n";}

	return($trim);		# いい感じに削った場合は１を返す。削らなかった場合は０を返す。
}

open REP, '>'.$report_file;
my $shortreads = 0;
my $trimbases = 0;
my $totalbases = 0;
for(my $f=0;$f<$files;$f++){			# メインループ。対象ファイルごとにやっていく
	print 'trimming '.$file[$f];
	open IN, $file[$f];	
	open OUT, '>'.$outfile[$f];
	my $fastq_count = 0;
	my $seq_count = 0;
	my $trimseqs = 0;
	$trimbases = 0;
	$totalbases = 0;
	$shortreads = 0;
	LOOP:while(my $name = <IN>){		# fastq を入力に取る。Nanopore のfastq は途中に改行は入ってないので４行読み込めばいい
		$name =~ s/\r|\n|\r\n//g;
		my $err = 0;
		if(substr($name,0,1) ne '@'){$err=1;}
		$seq = <IN>;
		$seq =~ s/\r|\n|\r\n//g;
		$totalbases += length($seq);
		my $plus = <IN>;
		$plus =~ s/\r|\n|\r\n//g;
		if($plus ne '+'){$err=1;}
		$qual = <IN>;
		$qual =~ s/\r|\n|\r\n//g;
		
		$trimseqs += edge();				# $seq, $qual の端を判定して処理する。関数にした
		if(length($seq)>$shortlimit){		# 短くなったやつは無かったことにする。もともと短かったやつもなかったことにする
			print OUT $name."\n";
			print OUT $seq."\n";
			print OUT $plus."\n";
			print OUT $qual."\n";
		}else{
			$shortreads++;
			$trimbases += length($seq);
		}
		
		if($err == 0){$seq_count++;}else{print REP 'wrong fastq format found in: '.$file[$f]."\n";last LOOP;}
		$fastq_count++;
		if($fastq_count>1000000){print REP 'too many reads in: '.$file[$f]."\n";last LOOP;}
		if($fastq_count % 1000 == 0){print '.';}
	}
	close OUT;
	close IN;
	print 'done'."\n";
	
	print REP $file[$f]."\n";	# レポート。
	my $shortreads_ratio = int(1000*$shortreads/$fastq_count)/10;	# 短すぎて消されたリードの数
	print REP 'length < '.$shortlimit.':'."\t".$shortreads_ratio.'% ('.$shortreads.'/'.$fastq_count.' reads)'."\n";
	my $trimseqs_ratio = int(1000*$trimseqs/$fastq_count)/10;		# トリムされたリードの数
	print REP 'trimed reads:'."\t".$trimseqs_ratio.'% ('.$trimseqs.'/'.$fastq_count.' reads)'."\n";
	my $final_seqs = $seq_count - $shortreads;						# 最終的に残ったリードの数
	my $error_ratio = int(1000*$final_seqs/$fastq_count)/10;
	print REP 'final reads:'."\t".$error_ratio.'% ('.$final_seqs.'/'.$fastq_count.' reads)'."\n";
	my $base_ratio = int(1000*$trimbases/$totalbases)/10;			# トリムされた塩基数。短いリードの削除も含む
	my $totalbasesm = int(($totalbases)/10000)/100;
	my $trimbasesm = int(($trimbases)/10)/100;
	print REP 'trimed bases:'."\t".$base_ratio.'% ('.$trimbasesm.' kb/'.$totalbasesm.' Mb)'."\n";
	my $finalbase_ratio = int(1000*($totalbases-$trimbases)/$totalbases)/10;			# トリムされた塩基数。短いリードの削除も含む
	my $finalbases = int(($totalbases-$trimbases)/10000)/100;
	print REP 'final bases:'."\t".$finalbase_ratio.'% ('.$finalbases.'/'.$totalbasesm.' Mb)'."\n";
	print REP "\n";

}
close REP;
close TRM;
