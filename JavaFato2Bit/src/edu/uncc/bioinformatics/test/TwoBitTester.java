package edu.uncc.bioinformatics.test;
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
import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import junit.framework.TestCase;

import org.junit.Test;

import edu.uncc.bioinformatics.DnaUtil;
import edu.uncc.bioinformatics.DuplicateSequenceException;
import edu.uncc.bioinformatics.FaFormatException;
import edu.uncc.bioinformatics.Fast2BitCompare;
import edu.uncc.bioinformatics.Fa2Bit;
import edu.uncc.bioinformatics.twobit.FileSizeExceedsTwoBitStandardException;
import edu.uncc.bioinformatics.twobit.TwoBitNameTooLongException;

public class TwoBitTester extends TestCase {
	/*
	@Test
	public void test() {
		fail("Not yet implemented");
	}
	*/
	@Test 
	public void testConverter(){
		String[] Input = {"TestData/chromosome1.fasta" , "TestData/chromosome2.fasta" , "TestData/chromosome3.fasta" };
		String Output = "TestData/java_fa_to_bit_output.2bit";
		
		long start = System.currentTimeMillis();
		try {
			
			Fa2Bit.Fastto2BitConvert( Input ,Output, false);
			
		} catch (FileSizeExceedsTwoBitStandardException e) {
			e.printStackTrace();
			fail("Exception");
		} catch (IOException e) {
			e.printStackTrace();
			fail("Exception");
		} catch (TwoBitNameTooLongException e) {
			e.printStackTrace();
			fail("Exception");
		} catch (FaFormatException e) {
			e.printStackTrace();
			fail("Exception");
		} catch (DuplicateSequenceException e) {
			e.printStackTrace();
			fail("Exception");
		}
		long stop = System.currentTimeMillis();
		System.out.println("Time: " + (stop - start) );
		
		String input = "";
		for( int i = 0; i < Input.length; i++){
			input += "file://" + new File(Input[i]).getAbsolutePath() + ( i == Input.length - 1 ? "" : "@" );
		}
		String [] args = {input, "file://" + new File(Output).getAbsolutePath() };
		Fast2BitCompare.main( args );
		assert( Fast2BitCompare.TestPassed );
	}
	@Test 
	public void testPackDna4(){
		DnaUtil.initNtChars();
		DnaUtil.initNtVal();
		{
			//test proper packing of bits into one byte
			char[] set = { 'a', 't', 'g', 'c' };
			Map<Character, String > map = new  HashMap<Character, String>();
			map.put( 'a', "10" );
			map.put( 't', "00" );
			map.put( 'c', "01" );
			map.put( 'g', "11" );
			
			for( int i = 0; i < 4; i++){
				for( int j = 0; j < 4; j++){
					for( int k = 0; k < 4; k++){
						for( int m = 0; m < 4; m++){
							char[] sample = { set[i], set[j], set[k], set[m] };
							char packed = DnaUtil.packDna4( sample, 0 );
							char check = (char) Integer.parseInt(map.get(set[i])+map.get(set[j])+map.get(set[k])+map.get(set[m]), 2); ; 
							//System.out.format("check %hx packed %hx", check, packed);
							assertEquals( packed , check );				
						}
					}
				}
			}	
		}
		{
			//test index shift
			char[] set = { 'a', 't', 'g', 'c' };
			Map<Character, String > map = new  HashMap<Character, String>();
			map.put( 'a', "10");
			map.put( 't', "00");
			map.put( 'c', "01");
			map.put( 'g', "11");
			
			for( int i = 0; i < 4; i++){
				for( int j = 0; j < 4; j++){
					for( int k = 0; k < 4; k++){
						for( int m = 0; m < 4; m++){
							char[] sample = { 'a', 't', 'a', set[i], set[j], set[k], set[m] };
							char packed = DnaUtil.packDna4(sample, 3);
							char check = (char) Integer.parseInt(map.get(set[i])+map.get(set[j])+map.get(set[k])+map.get(set[m]), 2); ; 
							//System.out.format("%hx", check);
							assert( packed == check );				
						}
					}
				}
			}
		}
		
	}

}
