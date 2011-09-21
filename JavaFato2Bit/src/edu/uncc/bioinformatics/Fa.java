package edu.uncc.bioinformatics;

import java.io.IOException;
import java.nio.CharBuffer;

public class Fa {

	
	/* Read in DNA or Peptide FA record in mixed case.   Allow any upper or lower case
	 * letter, or the dash character in.
	 * 
	 * boolean faMixedSpeedReadNext(struct lineFile *lf, DNA **retDna, int *retSize, char **retName)
	 * */
	public static boolean faMixedSpeedReadNext( LineFile f, DnaSeq seq ) throws IOException, FaFormatException{
		char c;
		int bufIx = 0;
		String name = null; //the name of teh fa sequence.
		int lineSize, i;
		String line;

		//this could be necesary 
		//dnaUtilOpen();

		/* Read first line, make sure it starts with '>', and read first word
		 * as name of sequence. */
		if ( (line = f.lineFileNext()) == null ){  //f, &line, &lineSize
		    //*retDna = NULL;
		    //*retSize = 0;
		    return false;
		}
		if (line.charAt(0) == '>'){
		    name = line.substring(1).replaceAll(" ", "");
		    if (line.equals("") ){
		        System.err.println("Expecting sequence name after '>' " + f.getFile().getName());
		    	throw new FaFormatException();
		    }
		}else{
		    System.err.println( "Expecting '>' line " + f.getFile().getName() );
		    throw new FaFormatException();
		}
		/* Read until next '>' */
		CharBuffer buff = CharBuffer.allocate( 10000000);
		for (;;){
		    if ( (line = f.lineFileNext()) == null ) //!lineFileNext(lf, &line, &lineSize)
		        break;
		    if (line.charAt(0) == '>'){
		    	f.lineFileReuse(line);
		    	break;
			}
		    for( char m : line.toCharArray()){
		    	if ( Character.isLetter(m) || m == '-'){
				    buff.put(m);
		    	}
		    }
		}
		//if (bufIx >= faFastBufSize)
		//    expandFaFastBuf(bufIx, 0);
		//faFastBuf[bufIx] = 0;
		//*retDna = faFastBuf;
		//*retSize = bufIx;
		char[] dna = new char[buff.position()];
		char[] c_buff = buff.array();
		for( int k = 0; k < dna.length; k++){
			dna[k] = c_buff[k];
		}
		seq.setName( name );
		seq.setDna( dna );
		
		//if (bufIx == 0){
		//    System.err.println("Invalid fasta format: sequence size == 0 for element " + name);
		//}
	
		return true;
	}

}
