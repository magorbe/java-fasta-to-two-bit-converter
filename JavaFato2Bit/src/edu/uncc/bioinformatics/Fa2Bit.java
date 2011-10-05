package edu.uncc.bioinformatics;
/** 
Copyright (C) 2011  Jeremy Villalobos

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

This program was a conversion from C languge.  The original code is from 
http://code.google.com/p/zinba/
**/
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
import java.io.BufferedOutputStream;

public class Fa2Bit {

	public static final String InputArg = "--input";
	public static final String OutputArg = "--output";
        public static final String MaskedArg = "--no-masked";
        public static final String ValidateArg = "--validate";
        public static final String HelpArg = "--help";
        public static final String LicenseArg = "--license";
        
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if( args.length < 2 ){
			System.out.println("Not enough arguments");
			printUsage();
			System.exit(1);
		}
		String[] Input = null;
		String Output = null;
                boolean Validate = false;
                //this false does translate soft and hard masked values.
                boolean NoMasked = false;
		for( int i = 0 ; i < args.length ; i++){
			if( args[i].equals(InputArg) ){
				String input = args[++i];
				Input = input.split("@");
			}else if (args[i].equals(OutputArg) ){
				Output = args[++i];
			}else if( args[i].equals( MaskedArg)){
                            NoMasked = true; //Boolean.parseBoolean( args[++i]);
                        }else if( args[i].equals(ValidateArg) ){
                            Validate = true;
                        }else if( args[i].equals( HelpArg)){
                            printUsage();
                            return;
                        }else if( args[i].equals( LicenseArg ) ) {
                            System.out.println(license);
                            return;
                        }
                        
		}
		if( Input == null ){
			System.out.println("Input is null.");
			printUsage();
			System.exit(2);
		}
                if( Output == null){
                    if( Input[0].endsWith(".fa")){
                        Output = Input[0].replace(".fa", ".2bit");
                    }else if(Input[0].endsWith( ".fasta" ) ) {
                        Output = Input[0].replace(".fa", ".2bit");
                    }
                    System.out.println("Output is null... setting the output file to " + Output );
                }
                if( NoMasked){
                    System.out.println("Converting masked nucleotide to proper letter");
                }else{
                    System.out.println("Converting masked nucleotide to n");
                }
		try {
			Fastto2BitConvert( Input , Output, NoMasked);
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
		if( Validate ){
	        File out = new File( Output );
	        Fast2BitCompare comp = new Fast2BitCompare();
	        boolean result = comp.validate( Input, "file://" + out.getAbsolutePath() );
	        if( !result ){
	            System.out.println("Test failed" );   
	        }
		}
	}
	
	public static void printUsage(){
		System.out.println("\n\tFato2Bit --input fastaone.fasta@fastwo.fasta@fastan.fasta --output out.2bit ");
		System.out.println("\t\t--input is required. At [@] is a special characer.  The fasta files should not include this character in their name.");
                System.out.println("\t\t--output can be used to specify the name of the 2bit output.  Otherwise the name of the first fa file is used.");
                System.out.println("\t\t--validate will check the fa and the 2bit files after creation.  The sequence and nucleotide position is printed if a check fails ");
                System.out.println("\t\t--no-masked : lower case nucleotide letters will be interpreted as n. instead of their respective nucleotide value.");
                System.out.println("\t\t\tNotice that the conversion is quiet fast without validation.");
                System.out.println("\t\t--help prints this message. \n");
                System.out.println("\t\t--license prints license. \n");
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
	public static void Fastto2BitConvert( String[] input, String output, boolean no_masked) throws FileSizeExceedsTwoBitStandardException, IOException
																			, TwoBitNameTooLongException, FaFormatException, DuplicateSequenceException{
		
		DnaUtil.initNtChars();
		DnaUtil.initNtVal();
		
		boolean noMaskFT = no_masked;
		boolean stripVersion = true;
		boolean ignoreDups = true;
		//create class called twoBit
		List<TwoBit> twoBitList = new ArrayList<TwoBit>();  //I will use a Java list instead of implementing a new link list.
		TwoBit twoBit; 
		Map<String, String> uniqHash = new HashMap<String, String>();
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
			System.out.println("Working on fa: " + fasta_file);
			Fa fa = new Fa();
			while( fa.faMixedSpeedReadNext( lf, seq)){ // (faMixedSpeedReadNext(lf, &seq.dna, &seq.size, &seq.name)
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
				/*This hash is only used to check for repetition of the sequence.  Therefore, I am only storing the 
				 * string twice instead of storeing the sequence, which produces a HUGE memory leak.*/
				uniqHash.put(seq.getName(), seq.getName() );
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
		//free up the big buffer memory.
		//Fa.invalidateBuffer();
		//I don't know why we need to revers the list.  but here we go...
		List<TwoBit> revers = new ArrayList<TwoBit>();
		for( int k = twoBitList.size() - 1; k >= 0; k--){
			revers.add( twoBitList.get( k ));
		}
		f = new File( output );
		
                BufferedOutputStream stream = new BufferedOutputStream( new DataOutputStream( new FileOutputStream( f ) ));
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
		

        static CharBuffer ans = CharBuffer.allocate(10000000);
	/* Convert non ACGT characters to N. */
	public static CharBuffer unknownToN(CharBuffer s){
		char c;
		ans.clear();
		char[] dna = s.toString().toCharArray();
		int i;
		for (i=0; i< dna.length; ++i){
		    c = dna[i];
                    char m = c;
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
                    char m = DnaUtil.ntChars[(int)buff[i]];
		    if ((c = m) == 0){
		    	c = 'n';
		    }
                 
		    buff[i] = c;
	    }
		
		return buff;
	}
	public static final String license = "Copyright (C) 2011  Jeremy Villalobos\n"+
            "Bioinformatics Department, University of North Carolina Charlotte.\n"+
            "This program is free software: you can redistribute it and/or modify\n"+
            "it under the terms of the GNU General Public License as published by\n"+
            "the Free Software Foundation, either version 3 of the License, or \n"+
            "(at your option) any later version.\n"+
            "\n"+
            "This program is distributed in the hope that it will be useful,\n"+
            "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"+
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"+
            "GNU General Public License for more details.\n"+
            "\n"+
            "You should have received a copy of the GNU General Public License\n"+
            "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"+
            "\n"+
            "This program was a conversion from C languge.  The original code is from \n"+
            "http://code.google.com/p/zinba/";
	
}
