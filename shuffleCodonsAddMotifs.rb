#!/usr/bin/ruby
require "DNATOOLBOX.rb"
require "HOME.rb"
inFP=ARGV[0];
outFP = ARGV[1];
motif = ARGV[2];
codonFP = ARGV[3];
if codonFP.nil?
	codonFP=HOME+"/genomes/sc/Genetic_Code.txt";
end
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

def fillInRest(aaSeq, setCodons, defaultStates, codonCount, aa2codon, codon2aa)
	0.upto(aaSeq.length-1){|i|
		if setCodons[i]==nil
			if codonCount[defaultStates[i]]>=1
				codonCount[defaultStates[i]]-=1;
				setCodons[i]=defaultStates[i];
			end
		end
	}
	#now that all defaults restored, fill in the others with the remaining codons.
	0.upto(aaSeq.length-1){|i|
		if setCodons[i]==nil
			foundRep=false;
			aa2codon[aaSeq[i]].each{|replacementCodon|
				if codonCount[replacementCodon]>=1
					foundRep=true;
					codonCount[replacementCodon]-=1;
					setCodons[i]=replacementCodon;
					break;
				end
			}
			if !foundRep
				raise("No more codons for: "+aaSeq[i]+" at "+i.to_s+"\n");
			end
		end
	}
	
	translated = translate(setCodons.join(""), codon2aa);
	if translated!=aaSeq.join("")
		raise("In and out AA seqs don't match: \n"+aaSeq.join("")+"\n"+translated+"\n");
	else
		#print("Same AA Seq!\n"+aaSeq.join("")+"\n"+translated+"\n");
	end
	return setCodons;
end

def getTheseShuffledOrfs(codon2aa,aa2codon,motif,curOrf,curOrfName,outFile)
	curOrfName=">"+curOrfName;
	if curOrf.length()%3!=0;
		raise("ORF NOT DIVISIBLE BY 3!!\n");
	end
	#print(curOrf+"\n");
	codonCount = Hash.new();
	aaSeq=[];
	setCodons=[];
	defaultStates=[];
	0.upto((curOrf.length()/3)-1){|i| #don't have to skip first and last since Met has only one codon and there is only one stop in an ORF
		curCodon = curOrf[i*3,3];
		defaultStates.push(curCodon);
		setCodons.push(nil);
		#print (curCodon+"\n");
		if codonCount[curCodon]==nil
			codonCount[curCodon]=0;
		end
		codonCount[curCodon]+=1;
		aaSeq.push(codon2aa[curCodon]);
	}
	codon2aa.keys.each{|codon|
		codonCount[codon]=0 if codonCount[codon]==nil;
	}
	# start on greedy algorithm.  see notes 2009/10/08-09
	#make DFA transitions
	offset=0; # what offset are we trying to start the motif at
	numIn=0; # how many codons in are we in the motif matching?
	functionalCodonCount=codonCount.clone(); #this is used to temporarily decrease the codon counts when trying greedily to assemble a motif
	curAA=0; # ORF codon index
	lastDoneAA=0;
	codonPossToSet=[];
	codonsToSet=[];
	foundOne=false
	oneIFound=nil;
	motifCount=0;
	
	while curAA<=aaSeq.length-1 #{ for each amino acid and each possible offset
		aa = aaSeq[curAA];
		if numIn==0
			fOffset=offset;
		else 
			fOffset=0;
		end
		mSt = [3*numIn-offset, 0].max;
		if fOffset==0
			mLen = [motif.length()-mSt,3].min();
			
		else
			mLen = [motif.length()-mSt,3-fOffset].min();
		end
		
		lookingFor = motif[mSt, mLen];
		#print("looking for: "+lookingFor+" "+curAA.to_s+" "+offset.to_s+" "+numIn.to_s+"\n");
		aa2codon[aa].each{|codon|
			if functionalCodonCount[codon]>=1 #still some of this codon remaining
				#print("  comparing to: "+codon[fOffset,mLen]+"...");
				if lookingFor==codon[fOffset,mLen]
					#print("good!\n");
					foundOne=true;
					oneIFound=codon;
					break;
				else
					#print("bad\n");
				end
			end
		}
		if foundOne #{
			#print("found one: "+oneIFound+" "+offset.to_s+" "+numIn.to_s+"\n");
			foundOne=false;
			functionalCodonCount[oneIFound]-=1;
			codonPossToSet.push(curAA);
			codonsToSet.push(oneIFound);
			curAA+=1;
			if mSt+mLen==motif.length #{
				#print("FULL HIT: ");
				#you got the full motif, finalize all in decrementList
				codonCount=functionalCodonCount;
				functionalCodonCount = codonCount.clone();
				0.upto(codonPossToSet.length-1){|i|
					setCodons[codonPossToSet[i]]=codonsToSet[i];
					if aaSeq[codonPossToSet[i]]!=codon2aa[codonsToSet[i]]
						raise("HALLELULEAH: "+codon2aa[codonsToSet[i]]+" "+aaSeq[codonPossToSet[i]]+" "+codonsToSet[i]+"\n");
					end
					#print(codonsToSet[i]);
				}
				#print("\n");
				codonPossToSet=[];
				codonsToSet=[];
				lastDoneAA=curAA;
				numIn=0;
				offset=0;
				motifCount+=1;
			else #}{
				numIn+=1;
			end #}
		else #}{   !foundOne
			if offset==2 #could not find any candidates for this aa, move on
				offset=0;
				lastDoneAA+=1;
			else # try to add the motif at an offset of two instead
				offset+=1;
			end
			#reset the motif matching
			curAA=lastDoneAA;
			numIn=0;
			codonPossToSet=[];
			codonsToSet=[];
			functionalCodonCount=codonCount.clone(); # reset to have the previous codon changes reversed
		end #}
	end #}
	#print("\n"+aaSeq.join()+"\n");
	#now fill in the others, with their default values if possible
	setCodons2=setCodons.clone();
	codonCount2=codonCount.clone();
	
	setCodons =  fillInRest(aaSeq, setCodons, defaultStates, codonCount, aa2codon, codon2aa);
	
	maximalSeq = setCodons.join("");
	actualCount = countMotif(maximalSeq, motif);
	if motifCount>actualCount
		raise("OH GOD! this should never happen!\n");
	elsif motifCount<actualCount
		print("Bonus!\n")
	end
	origSeq=defaultStates.join("");
	defaultCount=countMotif(origSeq, motif);
	outFile.print(curOrfName+" DEFAULT #"+motif+"="+defaultCount.to_s()+"\n"+origSeq+"\n");
	outFile.print(curOrfName+" MAX #"+motif+"="+actualCount.to_s()+"\n"+maximalSeq+"\n");
	lastCount = actualCount;
	i=0; #aaPosition
	j=0; #codon
	while i<=setCodons2.length-1
		if setCodons2[i]!=nil
			#print(i.to_s+" "+j.to_s+"\n");
			thisSetCodons = setCodons2.clone();
			thisCodonCount = codonCount2.clone();
			#change curCodon
			curCodon = setCodons2[i];
			if aa2codon[aaSeq[i]][j]==nil
				i+=1;
				j=0;
				next;
			elsif codonCount2[aa2codon[aaSeq[i]][j]]==0 || aa2codon[aaSeq[i]][j]==curCodon
				j+=1;
				next;
			end
			repCodon = aa2codon[aaSeq[i]][j];
			thisCodonCount[repCodon]-=1;
			thisCodonCount[curCodon]+=1;
			thisSetCodons[i]=repCodon;
			
			setCodons =  fillInRest(aaSeq, thisSetCodons, defaultStates, thisCodonCount, aa2codon, codon2aa);
			curSeq = setCodons.join("");
			newCount = countMotif(curSeq, motif);
			if newCount<lastCount
				lastCount = newCount;
				setCodons2[i]=repCodon;
				codonCount2[repCodon]-=1;;
				codonCount2[curCodon]+=1;;
				#print this one and save state
				#print("Found one that messes things up!\n");
				if newCount>defaultCount
					outFile.print(curOrfName+" MORE #"+motif+"="+newCount.to_s()+"\n"+curSeq+"\n");
				elsif newCount==defaultCount
					outFile.print(curOrfName+" EQUAL #"+motif+"="+newCount.to_s()+"\n"+curSeq+"\n");
				else #less
					outFile.print(curOrfName+" LESS #"+motif+"="+newCount.to_s()+"\n"+curSeq+"\n");
				end
				setCodons2=thisSetCodons;
				i+=1;
				j=0;
			else
				j+=1;
			end
		else
			i+=1;
		end
	end
	
end

outFile = File.open(outFP, "w");
curOrf=nil;
curOrfName="";
File.foreach(inFP){|line|
	line.chomp!();
	if line==nil||line==""
		next;
	elsif line[0,1]==">"
		if !curOrf.nil?
			print(curOrfName+"\n");
			getTheseShuffledOrfs(codon2aa,aa2codon,motif,curOrf,curOrfName,outFile)
		end
		curOrf="";
		curOrfName = line[1,line.length()-1];
	else
		curOrf+=line;
	end
}
if !curOrf.nil?
	print(curOrfName+"\n");
	getTheseShuffledOrfs(codon2aa,aa2codon,motif,curOrf,curOrfName,outFile)
end

outFile.close();
