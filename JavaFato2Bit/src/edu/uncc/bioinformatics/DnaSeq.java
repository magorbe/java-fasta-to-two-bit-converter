package edu.uncc.bioinformatics;

import java.nio.CharBuffer;
import java.util.ArrayList;

public class DnaSeq {
	//Use a Java List instead of implementing a linked list
	  //struct dnaSeq *next;  /* Next in list. */
	String Name;
	/**
	 * Stores the DNA data
	 * This should take care of buffer micromanaging while keeping memory footpring small 
	 */
	char[] Dna;
	//size not needed the Java array has size
	//skip Bits 
	//    char *name;           /* Name of sequence. */
	//    DNA *dna;             /* Sequence base by base. */
	//    int size;             /* Size of sequence. */
	//    Bits* mask;           /* Repeat mask (optional) */
	public String getName() {
		return Name;
	}
	public void setName(String name) {
		Name = name;
	}
	public char[] getDna() {
		return Dna;
	}
	public void setDna( char[] dna) {
		Dna = dna;
	}
	public int size(){
		return Dna.length;
	}
}
