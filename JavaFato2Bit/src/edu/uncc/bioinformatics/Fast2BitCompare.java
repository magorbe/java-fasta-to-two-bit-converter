/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.uncc.bioinformatics;
/** 
Copyright (C) 2011 Jeremy Villalobos
This module uses CPL licensed code.  It is a separate, independent module from the 
Fa22Bit module.  It is included for convenience.

**/
import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.SeqSymmetry;
import com.affymetrix.genometryImpl.span.SimpleSeqSpan;
import com.affymetrix.genometryImpl.symloader. Fasta;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URI;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jfvillal
 */
public class Fast2BitCompare {
	public static boolean TestPassed = true;
        
    static URI twobit;
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String[] files = args[0].split("@");
        Fast2BitCompare comp = new Fast2BitCompare();
        //comp.validate( files , args[1]);
        
        if( !comp.validate( files, args[1] )){
            System.out.println("Test failed" );   
        }
    }
    /**
     * 
     * @param input
     * @param output
     * @return returns true if both files are the same.
     */
    public boolean validate( String[] input, String output ){
    	twobit = URI.create(output);
        long start = System.currentTimeMillis();
        try {    
            int CPU = 1;
            Worker[] w = new Worker[CPU];
            for( int i = 0 ; i < CPU ; i++){
                w[i] = new Worker( twobit, queue );
                w[i].setName("worker(" + i + ")");
                w[i].start();
            }
            
            for( int i = 0; i < input.length; i++){
               queue.put( new Bucket(true) );
               doOneFile( input[i]);
            }
           
            Bucket b = new Bucket( true);
            b.setDone(true );
            queue.put( b );
            
            for( int i = 0 ; i < CPU ; i++){
                w[i].join();
            }
            
            
        } catch (Exception ex) {
            Logger.getLogger(Fast2BitCompare.class.getName()).log(Level.SEVERE, null, ex);
        }
        long stop = System.currentTimeMillis();
        //System.out.println("Time:" + (stop - start) );
        return TestPassed;
        
    }
    
    ArrayBlockingQueue<Bucket> queue = new ArrayBlockingQueue<Bucket>(1);
    
    public void doOneFile(String file) throws FileNotFoundException, IOException, FaFormatException, Exception{
        LineFile lf = new LineFile( file );
        DnaSeq seq = new DnaSeq(); //dnaSeq struct clone
        System.out.println("Working on fa: " + file);
        Fa fa = new Fa();
        while( fa.faMixedSpeedReadNext( lf, seq) ){ // 
            //gets the sequence for each file.

            char[] fast_region = seq.getDna();
            //char[] bit_region = bit_reg.toCharArray();
            Bucket b = new Bucket( fast_region,  seq.Name    ) ;
            
            queue.put( b );
        }
        lf.close();
    }
   
}