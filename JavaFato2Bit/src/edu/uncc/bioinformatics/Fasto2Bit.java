package edu.uncc.bioinformatics;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.CharBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uncc.bioinformatics.twobit.FileSizeExceedsTwoBitStandardException;
import edu.uncc.bioinformatics.twobit.TwoBit;
import edu.uncc.bioinformatics.twobit.TwoBitNameTooLongException;

public class Fasto2Bit {

	public static final String InputArg = "--input";
	public static final String OutputArg = "--output";
			
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if( args.length < 4 ){
			System.out.println("Not enough arguments");
			printUsage();
			System.exit(1);
		}
		String[] Input = null;
		String Output = null;
		for( int i = 0 ; i < args.length ; i++){
			if( args[i].equals(InputArg) ){
				String input = args[++i];
				Input = input.split("@");
			}else if (args[i].equals(OutputArg) ){
				Output = args[++i];
			}
		}
		if( Input == null || Output == null){
			System.out.println("Either Input or Output are null.");
			printUsage();
			System.exit(2);
		}
		
		try {
			Fastto2BitConvert( Input , Output);
		} catch (FileSizeExceedsTwoBitStandardException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (TwoBitNameTooLongException e) {
			e.printStackTrace();
		} catch (FaFormatException e) {
			e.printStackTrace();
		} catch (DuplicateSequenceException e) {
			e.printStackTrace();
		}
		
		System.out.println( "Done");
		

	}
	
	public static void printUsage(){
		System.out.println("Fato2Bit --input fastaone.fasta@fastwo.fasta@fastan.fasta --output out.2bit ");
		System.out.println("All the arguments are required. At [@] is a special characer.  The fasta files should not include this character in their name.");
	}
	/**
	 * Used for non-recoginized encodings.  
	 * TODO this method may be taken out because we expect the format to be known, and the format sould be fasta.
	 */
	/* Convert non ACGT characters to N. */
	static char[] unknownToN(char[] s ){
		char c;
		int i;
		for (i=0; i< s.length; ++i){
		    c = s[i];
		    if (DnaUtil.ntChars[(int)c] == 0){
				if (Character.isUpperCase(c)){
				    s[i] = 'N';
				}else{
				    s[i] = 'n';
				}
			}
		}	
		return s;
	}
	/** 
	 * Convert inFiles in fasta format to outfile in 2 bit 
	 * format. 
	 * @throws TwoBitNameTooLongException 
	 * @throws IOException 
	 * @throws FileSizeExceedsTwoBitStandardException 
	 * @throws DuplicateSequenceException 
	 * @throws FaFormatException 
	 **/	
	public static void Fastto2BitConvert( String[] input, String output) throws FileSizeExceedsTwoBitStandardException, IOException
																			, TwoBitNameTooLongException, FaFormatException, DuplicateSequenceException{
		
		DnaUtil.initNtChars();
		DnaUtil.initNtVal();
		
		boolean noMaskFT = false;
		boolean stripVersion = true;
		boolean ignoreDups = true;
		//create class called twoBit
		List<TwoBit> twoBitList = new ArrayList<TwoBit>();  //I will use a Java list instead of implementing a new link list.
		TwoBit twoBit; 
		Map<String, DnaSeq> uniqHash = new HashMap<String, DnaSeq>();
		//struct twoBit *twoBitList = NULL, *twoBit;
		/**
		 * I don't know what i is being used for
		 */
		int i;
		// Us a Java map for this, figure ou why is the hash needed.
		//struct hash *uniqHash = newHash(18);
		File f;
		//FILE *f;

		//int inFileCount = RinFileCount[0];
		/**
		 * inFiles are the input fasta files
		 * outFiles is the output twobit file
		 */
		//char *inFiles=RinFiles[0];
		//char *outFile = RoutFile[0];
		//char *delim = "@";  //this is already done in the main function
		//char *ptr= NULL;	
		for( String fasta_file : input ){
			LineFile lf = new LineFile( fasta_file); //lineFileOpen( fasta_file, true); 
													 //implement the lineFileOpen routine
													 //TODO implement the lineFile class
			DnaSeq seq = new DnaSeq(); //dnaSeq struct clone
			while( Fa.faMixedSpeedReadNext( lf, seq)){ // (faMixedSpeedReadNext(lf, &seq.dna, &seq.size, &seq.name)
				if (seq.size() == 0)
				{
					System.err.println("Skipping item "+ seq.Name + " which has no sequence.\n");
					continue;
				}
				
				/* strip off version number */
				if (stripVersion)
				{
					int val = seq.getName().indexOf(".");
					if( val != -1){
						seq.setName(seq.getName().substring(0, val) );
					}
				}
		
				if ( uniqHash.containsKey(seq.getName()) ){
					if (!ignoreDups){
						System.err.println("Duplicate sequence name " + seq.getName());
						throw new DuplicateSequenceException();
					}else{
						continue;
					}
				}
				uniqHash.put(seq.getName(), seq);
				//hashAdd(uniqHash, seq.name, NULL);
				if (noMaskFT){
					seq.setDna( faToDna( seq.getDna()) );
				}else{
					seq.setDna( unknownToN( seq.getDna() ) );
				}
				twoBit = TwoBit.twoBitFromDnaSeq(seq, !noMaskFT);
				twoBitList.add( twoBit );
				//slAddHead(&twoBitList, twoBit);
			}
			lf.close();
		}
		//I don't know why we need to revers the list.  but here we go...
		List<TwoBit> revers = new ArrayList<TwoBit>();
		for( int k = twoBitList.size() - 1; k >= 0; k--){
			revers.add( twoBitList.get( k ));
		}
		f = new File( output );
		DataOutputStream stream = new DataOutputStream( new FileOutputStream( f ) );
				//mustOpen(outFile, "wb"); //Open binary file
		TwoBit.twoBitWriteHeader(revers, stream);
		for( TwoBit two_bit : revers ){
			TwoBit.twoBitWriteOne(two_bit, stream);
		}
		/*for (twoBit = twoBitList; twoBit != NULL; twoBit = twoBit->next){
		    twoBitWriteOne(twoBit, f);
		}
		carefulClose(&f);*/
		stream.close();
		
	}
		

	/* Convert non ACGT characters to N. */
	public static CharBuffer unknownToN(CharBuffer s){
		char c;
		CharBuffer ans = CharBuffer.allocate(10000000);
		char[] dna = s.toString().toCharArray();
		int i;
		for (i=0; i< dna.length; ++i){
		    c = dna[i];
		    if (DnaUtil.ntChars[(int)c] == 0){
		    	if ( Character.isUpperCase(c)){
		    		ans.put('N');
		    	}else{
		    		ans.put('n');
		    	}
			}else{
				ans.put(c);
			}
	    }
		return ans;
	}
	
	/* Convert possibly mixed-case DNA to lower case.  Also turn
	 * any strange characters to 'n'.  Does not change size.
	 * of sequence. */
	public static char[] faToDna( char[] buff){
		int i;
		char c;
		for (i=0; i< buff.length; ++i){
		    if ((c = DnaUtil.ntChars[(int)buff[i]]) == 0){
		    	c = 'n';
		    }
		    buff[i] = c;
	    }
		
		return buff;
	}
	
	
}
