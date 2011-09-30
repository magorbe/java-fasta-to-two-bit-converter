package edu.uncc.bioinformatics.twobit;

import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;

import edu.uncc.bioinformatics.DnaSeq;
import edu.uncc.bioinformatics.DnaUtil;
import edu.uncc.bioinformatics.Sig;
import java.io.BufferedOutputStream;

/**
 * Two bit representation of DNA
 * @author jfvillal
 *
 */
public class TwoBit {
	//bits32 32 unsigned bits are replaced with an integer in Java
	TwoBit next;	/* Next sequence in list (may be taken out and use a List from outside the class*/
	String name;			/* Name of sequence. */
	char[] UbyteData;		/* DNA at two bits per base. */
	int size;		/* Size of this sequence. */
	int nBlockCount;		/* Count of blocks of Ns. */
	int[] nStarts;		/* Starts of blocks of Ns. */
	int[] nSizes;		/* Sizes of blocks of Ns. */
	int maskBlockCount;	/* Count of masked blocks. */
	int[] maskStarts;		/* Starts of masked regions. */
	int[] maskSizes;		/* Sizes of masked regions. */
	int reserved;		/* Reserved for future expansion. */
	
	/** 
	 * Return size when packed, rounding up. 
	 * 
	 **/
	static int packedSize(int unpackedSize){
		return ((unpackedSize + 3) >> 2);
	}
	
	
	private static int Segment;
	//*this code is not in use
	private static class worker extends Thread{
		char[] ubyte;
		char[] dna;
		int start;
		int end;
		public worker( char[] u, char[] c, int s, int e){
			ubyte = u;
			dna = c;
			start = s;
			end = e;
			//System.out.println( "start: " + start  + "div: " + start%4 + " stop: " + end  + " div: " + end%4);
		}
		public void run(){
			for( int i = start; i < end; i+=4){
				ubyte[i/4] = DnaUtil.packDna4(dna, i);
			}
		}
	}
	
	/* Convert dnaSeq representation in memory to twoBit representation.
	 * If doMask is true interpret lower-case letters as masked. */
	public static TwoBit twoBitFromDnaSeq( DnaSeq seq, boolean doMask){
		int ubyteSize = packedSize(seq.size());
		char[] ubyte_pt; 		//UBYTE *pt;
		
		ubyte_pt = new char[ubyteSize];
		
		TwoBit twoBit = new TwoBit(); // struct twoBit *twoBit;
		
		char dna_last4[] = new char[4];	/* Holds few bases. */
		char[] dna;  //		DNA *dna;
		int  end;
	
		/* Allocate structure and fill in name. */
		//AllocVar(twoBit);
		//pt = AllocArray(twoBit->data, ubyteSize);
		twoBit.setName( seq.getName()); //->name = cloneString(seq->name);
		twoBit.size = seq.size();
	
		/* Convert to 4-bases per byte representation. */
		//seq.getDna().capacity()-seq.getDna().remaining();
		//dna = seq.getDna().arrayOffset()Array.seq.getDna().array(); //string of dns chars
		dna = seq.getDna();
		
		end = dna.length - 4; //I don't know why -4
		
		int CPU = 4;
		int seg = end / CPU; // 4 cores
		
		/**
		 * packDna4 passes the array with a start index
		 * and the function should return a character, which we then 
		 * put into ubyte_pt.
		 * 
		 * this code did not improve anything.
		 */
		/*worker[] workers = new worker[CPU];
		int rem;
		for( int i = 0; i < CPU; i++){
			int start = i*seg - ((i*seg) % 4 );
			int stop = (i+1)*seg - (((i+1)*seg) % 4);
			
			if( i == CPU - 1){
				workers[i]  = new worker(ubyte_pt, dna, start , end );	
			}else{
				workers[i]  = new worker(ubyte_pt, dna, start, stop );
			}
			workers[i].start();
		}
		for( int i = 0; i < CPU; i++){
			try {
				workers[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}*/
		
		int i;
		for ( i = 0; i<end; i += 4){
			
		    ubyte_pt[i/4] = DnaUtil.packDna4(dna, i);
		    
		}
	
		/* Take care of conversion of last few bases. */
		dna_last4[0] = dna_last4[1] = dna_last4[2] = dna_last4[3] = 'T'; //00
		for( int j = i ; j < seq.size(); j++){
			dna_last4[j-i] = dna[j];
		}
		ubyte_pt[i/4] = DnaUtil.packDna4( dna_last4, 0 );
		
		//memcpy(last4, dna+i, seq->size-i);
		//*pt = packDna4(last4);
	
		/* The twobit format has indexes to quickly seek to areas of the dna.
		 * This part finds the blocks by looking for the pattern 'NN'
		 * Each block is then indexed into nStarts and its size is stored in 
		 * nSizes
		 * Deal with blocks of N. 
		 * */
		twoBit.countBlocksOfN( dna );
		if (twoBit.getNBlockCount() > 0){
			twoBit.nSizes = new int[ twoBit.getNBlockCount() ];
			twoBit.nStarts = new int[ twoBit.getNBlockCount() ];
			twoBit.storeBlocksOfN( dna );
		    //storeBlocksOfN(dna, seq->size, twoBit->nStarts, twoBit->nSizes);
		}
	
		/* Deal with masking */
		if (doMask){
			twoBit.countBlocksOfLower( dna );
		    //twoBit->maskBlockCount = countBlocksOfLower(dna, seq->size);
			if( twoBit.getMaskBlockCount() > 0){
				twoBit.maskStarts = new int[twoBit.getMaskBlockCount()];
				twoBit.maskSizes = new int[twoBit.getMaskBlockCount()];
				//AllocArray(twoBit->maskStarts, twoBit->maskBlockCount);
				//AllocArray(twoBit->maskSizes, twoBit->maskBlockCount);
				twoBit.storeBlocksOfLower(dna);
				//storeBlocksOfLower(dna, seq->size, 
				//	twoBit->maskStarts, twoBit->maskSizes);
			}
		}
		//set the dna on the 2bit data structure
		twoBit.UbyteData = ubyte_pt;
		
		return twoBit;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	
	/****
	 * 
	 *  Count number of blocks of N's (or n's) in s.
	 *  where s is the dna sequence, and N is the blocks that will be
	 *  indexed for fast access by a twobit parser.
	 *  @param s the a char[] string of dna letters.
	 *  
	 ****/
	void countBlocksOfN( char[] s ){
		boolean isN, lastIsN = false;
		char c;
		int blockCount = 0;
	
		for (int i=0; i< s.length; ++i){
		    c = s[i];
		    isN = (c == 'n' || c == 'N');
		    if (isN && !lastIsN)
			++blockCount;
		    lastIsN = isN;
		}
		nBlockCount = blockCount;
	}
	/**
	 * returns the nBlockCount.  To compute nBlockCount use
	 * {@link #countBlocksOfN(char[])}.
	 * @return
	 */
	public int getNBlockCount(){
		return nBlockCount;
	}
	
	/****
	 * 
	 * Store starts and sizes of blocks of N's.
	 * 
	 * One N will delimite the end and start of a new block
	 *  
	 ****/
	void storeBlocksOfN( char[] s ){
		boolean isN, lastIsN = false;
		int startN = 0;
		char c;
		//count moves the nBlock index
		int st_count = 0;
		int si_count = 0;
		int i;
		for (i=0; i < s.length; ++i){
		    c = s[i];
		    isN = (c == 'n' || c == 'N');
		    if (isN){
				if (!lastIsN){
					startN = i;
				}
			}else{
				if (lastIsN){
						nStarts[st_count++] = startN;
						nSizes[si_count++] = i - startN;
				}
			}
		    lastIsN = isN;
		}
		if (lastIsN){
				nStarts[st_count] = startN;
				nSizes[si_count] = i - startN;
		}
	}
	
	/* Count number of blocks of lower case letters. */
	void countBlocksOfLower( char[] s ){
		int i;
		boolean isLower, lastIsLower = false;
		int blockCount = 0;
	
		for (i=0; i< s.length; ++i){
		    isLower = Character.isLowerCase(s[i]);
		    if (isLower && !lastIsLower)
			++blockCount;
		    lastIsLower = isLower;
		}
		maskBlockCount = blockCount;
	}
	public int getMaskBlockCount(){
		return maskBlockCount;
	}
	
	/**
	 *  Store starts and sizes of blocks of lower case letters. 
	 **/
	public void storeBlocksOfLower(char[] s ){
		int i;
		boolean isLower, lastIsLower = false;
		int startLower = 0;
		int st_count = 0;
		int si_count = 0;
		for (i=0; i<size; ++i){
		    isLower = Character.isLowerCase(s[i]);
		    if (isLower){
		    	if (!lastIsLower){
		    		startLower = i;
		    	}
			}else{
				if (lastIsLower){
					maskStarts[st_count++] = startLower;
				    maskSizes[si_count++] = i - startLower;
				}
			}
		    lastIsLower = isLower;
	    }
		if (lastIsLower){
		    maskStarts[st_count] = startLower;
		    maskSizes[si_count] = i - startLower;
	    }
	}
	
	/**
	 *  Write out header portion of twoBit file, including initial
	 * index 
	 * @param stream the binary file to which we are going to output the data
	 * @throws FileSizeExceedsTwoBitStandardException 
	 * @throws IOException 
	 * @throws TwoBitNameTooLongException 
	 **/
	public static void twoBitWriteHeader(List<TwoBit> twoBitList, BufferedOutputStream stream ) throws FileSizeExceedsTwoBitStandardException
																			, IOException, TwoBitNameTooLongException{
		int sig = Sig.twoBitSig;
		int version = 0;
		int seqCount = twoBitList.size();
		int reserved = 0;
		int offset = 0;
		long counter = 0; 	/*this is 64 bits*/
							/* check for 32 bit overflow */
	
		/* Write out fixed parts of header. (write binary to file) */
		stream.write( intToByteArray(sig) );
		stream.write( intToByteArray(version));
		stream.write( intToByteArray(seqCount));
		stream.write( intToByteArray(reserved));
		//
		//writeOne(f, sig);
		//writeOne(f, version);
		//writeOne(f, seqCount);
		//writeOne(f, reserved);
	
		/* Figure out location of first byte past index.
		 * Each index entry contains 4 bytes of offset information
		 * and the name of the sequence, which is variable length. */
		//offset = sizeof(sig) + sizeof(version) + sizeof(seqCount) + sizeof(reserved);
		offset = 16; //offset after acouting for sig, version, seqCount and reserved.
		
		//Check the name length on the sequences complies with standards.
		for( TwoBit twoBit : twoBitList){
			int nameLen = twoBit.getName().length();
			if( nameLen > 255 ){
				System.err.println("Name " + twoBit.getName() + " too long");
				throw new TwoBitNameTooLongException();
			}
			offset += nameLen + 1 + 4; //4 being the sizeof(bit32)
		}
	
		/* Write out index. */
		for( TwoBit twoBit : twoBitList){
			int size = TwoBit.twoBitSizeInFile( twoBit );
			//stream.write( twoBit.getName().getBytes(Charset.defaultCharset()) );
			TwoBit.writeString( twoBit.getName() , stream);
			stream.write( intToByteArray(offset) );
			offset += size;
			counter += size;
			if( counter > Integer.MAX_VALUE ){
				 System.err.println(
						 	"Error in faToTwoBit, index overflow at " + twoBit.getName() 
						 	+ ". The 2bit format "
			                + "does not support indexes larger than "
						 	+ Integer.MAX_VALUE / 1000000000+"Gb, \n"
			                + "please split up into smaller files.\n" 
			                );
				 throw new FileSizeExceedsTwoBitStandardException();
			}
		}
	}
	
	/****
	 * 
	 * Write a 255 or less character string to a file.
	 * This will write the length of the string in the first
	 * byte then the string itself.
	 * @throws TwoBitNameTooLongException 
	 * @throws IOException 
	 *  
	 ****/
	public static void writeString( String s, BufferedOutputStream stream) throws TwoBitNameTooLongException, IOException{
		char[] in = new char[s.length()];
		for( int i = 0; i < s.length(); i++){
			in[i] = s.charAt( i );
		}
		
		int len = in.length;
	
		if (len > 255){
		    	System.err.println("String too long in writeString (" + len + " chars):\n" + s);
		    	throw new TwoBitNameTooLongException();
		    	//len = 255;
		}
		
		//writeOne(f, bLen);
		stream.write(len);
		for( int i = 0 ; i < in.length; i++){
			stream.write( in[i] );
		}
		//mustWrite(f, s, len);
	}
	
	
	/**
	 * 
	 *  Figure out size structure will take in file. 
	 **/
	public static int twoBitSizeInFile(TwoBit twoBit){
		//for size, nBlockCount, maskBlockCount, and reserved we have 4 ints
		//so 4 * 4 = 16
		//nStarts, nSizes, maskStarts and maskSizes each would be 4 * length
		return packedSize( twoBit.size ) + 16 
				+ (twoBit.nStarts != null ? twoBit.nStarts.length    * 4 : 0) 
				+ (twoBit.nSizes != null ? twoBit.nSizes.length     * 4 : 0)
				+ (twoBit.maskStarts != null ? twoBit.maskStarts.length * 4 : 0 )
				+ (twoBit.maskSizes != null ? twoBit.maskSizes.length  * 4 : 0 );
		
		//this is the original c code
		/*return packedSize(twoBit->size) 
		+ sizeof(twoBit->size)
		+ sizeof(twoBit->nBlockCount)
		+ sizeof(twoBit->nStarts[0]) * twoBit->nBlockCount
		+ sizeof(twoBit->nSizes[0]) * twoBit->nBlockCount
		+ sizeof(twoBit->maskBlockCount)
		+ sizeof(twoBit->maskStarts[0]) * twoBit->maskBlockCount
		+ sizeof(twoBit->maskSizes[0]) * twoBit->maskBlockCount
		+ sizeof(twoBit->reserved); */
	}
	
	
	public static final byte[] intToByteArray(int value) {
        return new byte[] {
                (byte)(value >>> 24),
                (byte)(value >>> 16),
                (byte)(value >>> 8),
                (byte)value};
	}
	
	/** 
	 * 
	 * Write out one twoBit sequence to binary file. 
	 * Note this does not include the name, which is
	 * stored only in index. 
	 * @throws IOException 
	 *
	 **/
	public static void twoBitWriteOne(TwoBit twoBit, BufferedOutputStream f) throws IOException{
		f.write( intToByteArray( twoBit.size ) );
		//write the count for the n block count.
		f.write( intToByteArray( twoBit.nBlockCount ));
		//write all the sequence data.
		if( twoBit.nBlockCount > 0){
			for( int i = 0; i < twoBit.nBlockCount; i++){
				f.write( intToByteArray( twoBit.nStarts[i] ));
			}
			for( int i = 0; i < twoBit.nBlockCount; i++){
				f.write( intToByteArray( twoBit.nSizes[i]));
			}
		}
		//write the count for the mask block
		f.write( intToByteArray( twoBit.maskBlockCount ));
		//write all the sequence data.
		if( twoBit.maskBlockCount > 0 ){
			for( int i = 0; i < twoBit.maskBlockCount; i++){
				f.write( intToByteArray(  twoBit.maskStarts[i] ) );
			}
			for( int i = 0; i < twoBit.maskBlockCount; i++){
				f.write( intToByteArray( twoBit.maskSizes[i] ));
			}
		}
		f.write( TwoBit.intToByteArray(twoBit.reserved) );
		
		//must write ?
		for( int i = 0; i < twoBit.UbyteData.length; i++){
			//each char has one byte of inportant data.
			//with two bits for each of the nucleotides
			//because the char is 2 bytes, but we are not using
			//the upper byte, this method should produce the intended
			//result (write one bye for each char)
			f.write( twoBit.UbyteData[i] );
		}
		
		//below is some of the original c code.
		//writeOne(f, twoBit->size);
		//writeOne(f, twoBit->nBlockCount);
		/*if (twoBit->nBlockCount > 0){
		    fwrite(twoBit->nStarts, sizeof(twoBit->nStarts[0]), 
		    	twoBit->nBlockCount, f);
		    fwrite(twoBit->nSizes, sizeof(twoBit->nSizes[0]), 
		    	twoBit->nBlockCount, f);
		}*/
		//
		
		//writeOne(f, twoBit->maskBlockCount);
		/*if (twoBit->maskBlockCount > 0)
		    {
		    fwrite(twoBit->maskStarts, sizeof(twoBit->maskStarts[0]), 
		    	twoBit->maskBlockCount, f);
		    fwrite(twoBit->maskSizes, sizeof(twoBit->maskSizes[0]), 
		    	twoBit->maskBlockCount, f);
		    }
		    */
		
		//writeOne(f, twoBit->reserved);
		//mustWrite(f, twoBit->data, packedSize(twoBit->size));
	}
	
}
