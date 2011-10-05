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
