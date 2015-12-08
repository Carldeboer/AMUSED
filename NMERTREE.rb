class NMerTreeHash
	def initialize(s,a)
		@maxTreeSize=s;
		@alphabet=a
		@seqHashes=Hash.new();
		1.upto(s){|i|
			enumerateKMers(@alphabet, i).each{|mer| @seqHashes[mer]=0};
		}
	end
	def traverseSeqs(seqFile)
		seqTotLen=0;
		seqCount  =0;
		prevSeq = "";
		while curSeq = seqFile.nextPart()
			if curSeq=="\n"#done with the last seq
				traverse(prevSeq);
				seqCount+=1;
				prevSeq ="";
			else
				curSeqs = curSeq.split(/[^ATGC]+/);
				while curSeqs.length()>1
					curSeq = curSeqs.shift();
					seqTotLen+=curSeq.length();
					traverse(prevSeq + curSeq);
					prevSeq = "";
					# because there is a break, count this as a new sequence.
					seqCount+=1;
				end
				curSeq = curSeqs.shift();
				if curSeq.nil?
					traverse(prevSeq) if !prevSeq.nil?
					prevSeq = "";
				else
					seqTotLen+=curSeq.length();
					curSeq = prevSeq+curSeq;
					partialTraverse(curSeq)
					if curSeq.length >= @maxTreeSize-1
						prevSeq = curSeq[curSeq.length()-@maxTreeSize+1,@maxTreeSize-1];
					else
						prevSeq = curSeq;
					end
				end
			end
		end
		if prevSeq!=""
			traverse(prevSeq);
			seqCount+=1;
		end
		return [seqCount, seqTotLen];
	end
	def partialTraverse(seq)
		#print(seq)
		1.upto(@maxTreeSize){|i| # for each k mer
			0.upto(seq.length()-@maxTreeSize){|j| # for each starting position
				@seqHashes[seq[j,i]]+=1;
			}
		}
	end
	def traverse(seq)
		1.upto(@maxTreeSize){|i| # for each k mer
			0.upto(seq.length()-i){|j| # for each starting position
				@seqHashes[seq[j,i]]+=1;
			}
		}
	end
	def getAllNMers(l,alphabet)
		#return Hash[enumerateAllKMers(alphabet,l).map{|key| [key,@seqHashes[key]]}];
		return enumerateAllKMers(alphabet,l).map{|key| [key,@seqHashes[key]]};
	end
	def keys(i)
		retMe = [];
		@seqHashes.keys.each{|k| 
			if k.length==i
				retMe.push(k);
			end
		}
		return retMe
	end
	def [](k)
		return @seqHashes[k];
	end
	def []=(k,v)
		return @seqHashes[k]=v;
	end
end

class KMerDist
	def initialize(alphabet,theLen,maxTreeSize)
		@alphabet = alphabet;
		@theLen=theLen;
		@maxTreeSize=maxTreeSize;
		@seqHashes=Hash.new();
		enumerateAllKMers(@alphabet, @maxTreeSize).each{|mer| @seqHashes[mer]=[0]*@theLen};
	end
	def keys()
		return @seqHashes.keys();
	end
	def [](key)
		return @seqHashes[key];
	end
	def traverse(seq)
		1.upto(@maxTreeSize){|i| # for each k mer
			0.upto(@theLen-i){|j| # for each starting position
				@seqHashes[seq[j,i]][j]+=1;
			}
		}
	end
end
def enumerateKMers(alphabet,k)
	return [""] if k==0;
	allKMers=[]
	kMers = enumerateKMers(alphabet,(k-1));
	kMers.each(){|mer|
		alphabet.each(){|b|
			allKMers.push(mer+b);
		}
	}
	return allKMers;
end

def enumerateAllKMers(alphabet,k)
	allKMers=[];
	1.upto(k){|i|
		allKMers+=enumerateKMers(alphabet,i);
	}
	return allKMers;
end
def revcomp(sequence)
	return sequence.reverse.tr("ATGCMRWSYKVHDB", "TACGKYSWRMBDHV"); #added for bad genomes 20110128
end

class SeqReader
	def initialize(file, isFasta, toUpper)
		@file=file;
		@isFasta = isFasta;
		@toUpper = toUpper;
		if @isFasta
			line = @file.gets();
			if line[0,1]==">"
				@hasNext = true;
			else
				raise("File not in FASTA format!");
			end
		end
	end
	def processSeq(seq)
		if @toUpper
			return seq.upcase();
		else
			return seq;
		end
	end
	def nextFA()
		curSeq = ""
		while line = @file.gets()
			if line[0,1]==">"
				return processSeq(curSeq);
			else
				curSeq+=line.chomp();
			end
		end
		@hasNext=false;
		return processSeq(curSeq);
	end
	#return the next sequence part (e.g. line of a fasta file), or "\n" if at the end of the sequence, or nil if at EOF
	def nextPart()
		if @isFasta
			line= @file.gets();
			if line.nil?
				return nil;
			elsif line[0,1]==">"
				@hasNext=true;
				return "\n";
			else
				return processSeq(line.chomp());
			end
		else
			if @halfPart
				@halfPart=false;
				return "\n";
			else
				line = @file.gets();
				@halfPart=true;
				return nil if line.nil?; #EOF
				return processSeq(line.chomp());
			end
		end
	end
	def next()
		if @isFasta
			if @hasNext
				return self.nextFA();
			else
				return nil;
			end
		else
			line = @file.gets();
			return nil if line.nil?;
			return processSeq(line.chomp());
		end
	end
end


# the following are vestigial.  Using a hash turned out to beslightly more efficient.

class NMerTree
	def initialize(maxLen)
		@maxDepth=maxLen;
		@rootNode = NMTNode.new();
	end
	
	def traverse(string)
		#string = string+"$";
		0.upto(string.length-1){|i|
			cur = string[i,(string.length-i)];
			#print(cur+"\n")
			@rootNode.traverse(cur, @maxDepth);
		}
	end
	def getThisMer(mer)
		return @rootNode.getThisMer(mer);
	end
	def getRandomMer(n)
		allNmers = @rootNode.getNMers(n,"");
		outStr=""; #TODO: not finished, was supposed to give a random string of length n that matched the base content of the tree
		
	end
	def getNMers(n)
		return @rootNode.getNMers(n, "");
	end
	
	def getAllNMers(n, alphabet)
		return @rootNode.getTheseNMers(n, "", alphabet);
	end
	def getSummary(depth)
		return @rootNode.getSummary(depth);
	end
	def addGaps(numThreads) 
		if numThreads==1
			@rootNode.addGapsNT(@maxDepth) 
		else
			@rootNode.addGaps(@maxDepth,numThreads) 
		end
	end
	def test()
		@rootNode.test();
	end
	def toString()
		@rootNode.toString();
	end
	def clear()
		@rootNode = NMTNode.new();
	end
end

class NMTNode
	attr_accessor :hits, :children;
	def initialize()
		@hits = 0;
		@children = Hash.new();
	end
	def test()
		print(@children.keys.length.to_s+"\n");
		
	end
	def toString()
			curS = "";
			@children.keys.each(){|k|
				curS += k+":"+@children[k].hits.to_s()+"{"+@children[k].toString()+"},";
			}
			return curS;
	end
	def addSubtree(subtree)#add the contents of this subtree to myself
		@hits+=subtree.hits;
		subtree.children.keys.each{|b|
			@children[b]=NMTNode.new() if !@children.has_key?(b);
			@children[b].addSubtree(subtree.children[b]);
		}
	end
	def addGapsNT(remainingDepth)
		if remainingDepth>0
			gapNode=NMTNode.new();
			@children.delete("-");
			#now add children to the gap node by traversing the others and summing to get the gapped
			@children.keys.each(){|k|
				gapNode.addSubtree(@children[k]);
			}
			@children["-"]=gapNode;
			#traverse children too
			@children.keys.each(){|k|
				@children[k].addGapsNT(remainingDepth-1);
			}
			#puts @children.keys()
		end
	end
	def addGaps(remainingDepth, numThreads)
		if remainingDepth>0
			gapNode=NMTNode.new();
			@children.delete("-");
			#now add children to the gap node by traversing the others and summing to get the gapped
			threads =[];
			@children.keys.each(){|k|
				if numThreads>Thread.list.count {|t| t.status == "run"} #there are enough to start a thread for the children
					threads.push(Thread.new{
						gapNode.addSubtree(@children[k]);
					});
				else
					gapNode.addSubtree(@children[k]);
				end
			}
			threads.each(){|t| t.join()};
			@children["-"]=gapNode;
			#traverse children too
			threads =[];
			@children.keys.each(){|k|
				if numThreads>Thread.list.count {|t| t.status == "run"} #there are enough to start a thread for the children
					threads.push(Thread.new{
						@children[k].addGaps(remainingDepth-1,numThreads)
					});
				else
					@children[k].addGapsNT(remainingDepth-1);# do it myself.  If new CPUs become free during this time, the next pass will get it, or the other addGaps calls on other threads will fill them
				end
			}
			threads.each(){|t| t.join()};
			#puts @children.keys()
		end
	end
	def getThisMer(string)
		return @hits if string.length()==0;
		thisTrans=string[0,1];
		rest = string[1,(string.length-1)];
		return 0 if @children[thisTrans]==nil;
		return @children[thisTrans].getThisMer(rest);
	end
	def traverse(string, remainingDepth)
		@hits+=1;
		if remainingDepth==0 || string==""
			return
		end
		remainingDepth-=1;
		thisTrans=string[0,1];
		rest = string[1,(string.length-1)];
		if @children[thisTrans]==nil
			@children[thisTrans]=NMTNode.new()
		end
		@children[thisTrans].traverse(rest, remainingDepth);
	end
	
	def getNMers(n,stringSoFar)
		if n==0;
			#print("Terminal: "+stringSoFar+"\t"+@hits.to_s+"\n");
			return [[stringSoFar, @hits]];
		end
		n-=1;
		results=[];
		
		@children.keys.each{|char|
			#print("Recursing for "+stringSoFar+"+"+char+" "+n.to_s+"\n");
			thisRes = @children[char].getNMers(n, stringSoFar+char);
			thisRes.each{|cur|
				results.push(cur);
			}
		}
		return results;
	end

	def getSummary(depth)
		depth-=1;
		#if @children
		if depth==1
			#nonRecurse
		else
			#recurse
		end
	end

	def getTheseNMers(n,stringSoFar, alphabet)
		if n==0;
			#print("Terminal: "+stringSoFar+"\t"+@hits.to_s+"\n");
			return [[stringSoFar, @hits]];
		end
		n-=1;
		results=[];
		
		0.upto(alphabet.length-1){|i|
			char = alphabet[i];
			if @children[char]==nil
				#the rest of these are 0
				curBin =[0]*n;
				done =false;
				while !done
					curStr = char;
					curBin.each{|b|
						curStr = curStr+alphabet[b];
					}
					results.push([stringSoFar+curStr, 0]);
					#increment
					binPos = n-1;
					carry=1;
					while binPos>=0 && carry>0
						if curBin[binPos]==alphabet.length-1
							curBin[binPos]=0;
							binPos-=1;
						else
							curBin[binPos]+=1;
							carry-=1;
						end
					end
					if carry==1
						done=true;
					end
				end
				next;
			end
			#print("Recursing for "+stringSoFar+"+"+char+" "+n.to_s+"\n");
			thisRes = @children[char].getTheseNMers(n, stringSoFar+char, alphabet);
			thisRes.each{|cur|
				results.push(cur);
			}
		}
		return results;
	end
end
