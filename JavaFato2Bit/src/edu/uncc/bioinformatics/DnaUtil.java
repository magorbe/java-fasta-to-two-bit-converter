package edu.uncc.bioinformatics;
/**
 * Some code from dnautil may be placed here if necesary
 * @author jfvillal
 *
 */
public class DnaUtil {
	/* Numerical values for bases. */
	static final int T_BASE_VAL = 0;
	static final int  U_BASE_VAL = 0;
	static final int  C_BASE_VAL = 1;
	static final int  A_BASE_VAL = 2;
	static final int  G_BASE_VAL = 3;
	static final int  N_BASE_VAL = 4;
	/* Used in 1/2 byte representation. */
	static final int MASKED_BASE_BIT = 8;
	
	static int[] ntVal = new int[256];
	static int[] ntValLower = new int[256];	/* NT values only for lower case. */
	static int[] ntValUpper = new int [256];	/* NT values only for upper case. */
	static int[] ntVal5 = new int[256];
	static int[] ntValNoN = new int[256]; /* Like ntVal, but with T_BASE_VAL in place of -1 for nonexistent ones. */
	
	static char[] valToNt = new char[(N_BASE_VAL|MASKED_BASE_BIT)+1];

	/* convert tables for bit-4 indicating masked */
	static int[] ntValMasked = new int[256];
	static char[] valToNtMasked = new char[256];

	static boolean inittedNtVal = false;

	static void initNtVal()
	{
	if (!inittedNtVal)
	    {
	    int i;
	    for (i=0; i< ntVal.length; i++)
	        {
		ntValUpper[i] = ntValLower[i] = ntVal[i] = -1;
	        ntValNoN[i] = T_BASE_VAL;
		if (Character.isSpaceChar((char)i) || Character.isDigit((char)i) ){
		    ntVal5[i] = ntValMasked[i] = -1;
		}else{
		    ntVal5[i] = N_BASE_VAL;
		    ntValMasked[i] = (Character.isLowerCase((char)i) ? (N_BASE_VAL|MASKED_BASE_BIT) : N_BASE_VAL);
	            }
	    }
	    ntVal5['t'] = ntVal5['T'] = ntValNoN['t'] = ntValNoN['T'] = ntVal['t'] = ntVal['T'] = 
	    	ntValLower['t'] = ntValUpper['T'] = T_BASE_VAL;
	    ntVal5['u'] = ntVal5['U'] = ntValNoN['u'] = ntValNoN['U'] = ntVal['u'] = ntVal['U'] = 
	    	ntValLower['u'] = ntValUpper['U'] = U_BASE_VAL;
	    ntVal5['c'] = ntVal5['C'] = ntValNoN['c'] = ntValNoN['C'] = ntVal['c'] = ntVal['C'] = 
	    	ntValLower['c'] = ntValUpper['C'] = C_BASE_VAL;
	    ntVal5['a'] = ntVal5['A'] = ntValNoN['a'] = ntValNoN['A'] = ntVal['a'] = ntVal['A'] = 
	    	ntValLower['a'] = ntValUpper['A'] = A_BASE_VAL;
	    ntVal5['g'] = ntVal5['G'] = ntValNoN['g'] = ntValNoN['G'] = ntVal['g'] = ntVal['G'] = 
	    	ntValLower['g'] = ntValUpper['G'] = G_BASE_VAL;

	    valToNt[T_BASE_VAL] = valToNt[T_BASE_VAL|MASKED_BASE_BIT] = 't';
	    valToNt[C_BASE_VAL] = valToNt[C_BASE_VAL|MASKED_BASE_BIT] = 'c';
	    valToNt[A_BASE_VAL] = valToNt[A_BASE_VAL|MASKED_BASE_BIT] = 'a';
	    valToNt[G_BASE_VAL] = valToNt[G_BASE_VAL|MASKED_BASE_BIT] = 'g';
	    valToNt[N_BASE_VAL] = valToNt[N_BASE_VAL|MASKED_BASE_BIT] = 'n';

	    /* masked values */
	    ntValMasked['T'] = T_BASE_VAL;
	    ntValMasked['U'] = U_BASE_VAL;
	    ntValMasked['C'] = C_BASE_VAL;
	    ntValMasked['A'] = A_BASE_VAL;
	    ntValMasked['G'] = G_BASE_VAL;

	    ntValMasked['t'] = T_BASE_VAL|MASKED_BASE_BIT;
	    ntValMasked['u'] = U_BASE_VAL|MASKED_BASE_BIT;
	    ntValMasked['c'] = C_BASE_VAL|MASKED_BASE_BIT;
	    ntValMasked['a'] = A_BASE_VAL|MASKED_BASE_BIT;
	    ntValMasked['g'] = G_BASE_VAL|MASKED_BASE_BIT;

	    valToNtMasked[T_BASE_VAL] = 'T';
	    valToNtMasked[C_BASE_VAL] = 'C';
	    valToNtMasked[A_BASE_VAL] = 'A';
	    valToNtMasked[G_BASE_VAL] = 'G';
	    valToNtMasked[N_BASE_VAL] = 'N';

	    valToNtMasked[T_BASE_VAL|MASKED_BASE_BIT] = 't';
	    valToNtMasked[C_BASE_VAL|MASKED_BASE_BIT] = 'c';
	    valToNtMasked[A_BASE_VAL|MASKED_BASE_BIT] = 'a';
	    valToNtMasked[G_BASE_VAL|MASKED_BASE_BIT] = 'g';
	    valToNtMasked[N_BASE_VAL|MASKED_BASE_BIT] = 'n';

	    inittedNtVal = true;
	    }
	}
	/**
	 *  
	 *  packs four bases (2 bits per base) into one byte 
	 *  and returns it as a char
	 *  the Java char is 2 bytes.
	 *  
	 *  Pack 4 bases into a UBYTE
	 *  @param in the array of dna letters
	 *  @offset where to start to collect letters. 
	 **/
	public static char packDna4(char[] in, int offset){
		char out = 0;
		int bVal;
		for( int i = offset; i < (offset + 4); i++){
		    bVal = ntValNoN[(int)in[i]];
		    out <<= 2;
		    out += bVal;
		}
		
		return out;
	}
	
	
	/* A little array to help us decide if a character is a 
	 * nucleotide, and if so convert it to lower case. */
	static char[] ntChars = new char[256];
	static boolean initted = false;
	static void initNtChars(){
		if (!initted){
		    ntChars['a'] = ntChars['A'] = 'a';
		    ntChars['c'] = ntChars['C'] = 'c';
		    ntChars['g'] = ntChars['G'] = 'g';
		    ntChars['t'] = ntChars['T'] = 't';
		    ntChars['n'] = ntChars['N'] = 'n';
		    ntChars['u'] = ntChars['U'] = 'u';
		    ntChars['-'] = 'n';
		    initted = true;
	    }
	}
}

