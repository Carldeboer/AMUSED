#!/usr/bin/ruby
require "DNATOOLBOX.rb"
require "HOME.rb"
inFP=ARGV[0];#fasta file
outFP = ARGV[1];
codonFP = HOME+"/genomes/sc/Genetic_Code.txt";

aa2codon = Hash.new();
codon2aa = Hash.new();
File.foreach(codonFP){|line|
	line.chomp!();
	if line==nil||line==""
		next;
	end
	codon, aa = line.split("\t");

	if codon2aa[codon]==nil
		codon2aa[codon] = aa;
	else
		raise("codon maps to more than one place: "+codon+"\n");
	end

	if aa2codon[aa]==nil
		aa2codon[aa]=[];
	end
	aa2codon[aa].push(codon);
}

allORFs = [];
allLabels = [];
theOrf=""
lastLabel="";
maxLen = 0;
File.foreach(inFP){|line|
	line.chomp!();
	if line.nil?||line==""
		next;
	elsif line[0,1]==">"
		if theOrf.length>0
			allORFs.push(theOrf);
			maxLen = [maxLen,theOrf.length()].max();
			allLabels.push(lastLabel);
		end
		theOrf="";
		lastLabel=line[1,line.length()-1];
		next;
	end
	theOrf+=line;
}
allORFs.push(theOrf);
allLabels.push(lastLabel);



outFile = File.open(outFP, "w");
0.upto(allORFs.length()-1){|i|
	codons=Hash.new();
	codingSeq=[];
	0.upto(allORFs[i].length()/3-1){|k|
		curCodon=allORFs[i][k*3,3];
		raise("Something wrong: #{curCodon}, k=#{(k*3).to_s()}, orf len=#{allORFs[i].length()}") if curCodon.length()!=3;
		curAA=codon2aa[curCodon];
		codingSeq.push(curAA);
		codons[curAA]=[] if codons[curAA].nil?;
		codons[curAA].push(curCodon);
	}
	#now shuffle
	curSeq="";
	0.upto(codingSeq.length()-1){|k|
		curAA=codingSeq[k];
		thisOne = (rand()*codons[curAA].length()).floor();
		curSeq+=codons[curAA].delete_at(thisOne);
	}
	outFile.print(">"+allLabels[i]+"\n"+curSeq+"\n");
}
outFile.close();
