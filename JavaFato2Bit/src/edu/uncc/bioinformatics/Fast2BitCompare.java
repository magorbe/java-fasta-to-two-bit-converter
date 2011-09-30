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
        static HashMap<String,BioSeq> map;
        static TwoBit bit;
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String[] files = args[0].split("@");
        validate( files , args[1]);
        
        if( !validate( files, args[1] )){
            System.out.println("Test failed" );   
        }
    }
    /**
     * 
     * @param input
     * @param output
     * @return returns true if both files are the same.
     */
    public static boolean validate( String[] input, String output ){
        long start = System.currentTimeMillis();
        try {
            
            map = new HashMap<String, BioSeq>();
            
            
            twobit = URI.create( output );
            AnnotatedSeqGroup grp_bit = new AnnotatedSeqGroup("bgroup");
            //this twobit is used to load the sequences into the map.
            TwoBit bit = new TwoBit( twobit, "", grp_bit);
            
            
            List<BioSeq> bit_lst = bit.getChromosomeList();
            System.out.println( "List size: " + bit_lst.size() );
            
            //int div = bit_lst.size() / 4;
            //for( int )
            //bit_lst.subList( , );
            
            for( BioSeq seq : bit_lst ){
                map.put(seq.getID(), seq);    
            }
            int CPU = 2;
            worker[] w = new worker[CPU];
            for( int i = 0 ; i < CPU ; i++){
                w[i] = new worker();
                w[i].setName("worker(" + i + ")");
                w[i].start();
            }
            
            for( int i = 0; i < input.length; i++){
               doOneFile( input[i]);
            }
            
            done = true;
            
            
            synchronized( semaphore){
                
                semaphore.notifyAll();
            }
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
        
    static public void doOneFile(String file) throws FileNotFoundException, IOException, FaFormatException, Exception{
        LineFile lf = new LineFile( file );
        DnaSeq seq = new DnaSeq(); //dnaSeq struct clone
        System.out.println("Working on fa: " + file);
        
        while( Fa.faMixedSpeedReadNext( lf, seq)){ // 
            //gets the sequence for each file.

            

            char[] fast_region = seq.getDna();
            //char[] bit_region = bit_reg.toCharArray();
            bucket b = new bucket( fast_region,  seq.Name    ) ;
            synchronized( queue){
                queue.add( b );
            }
            synchronized( semaphore){
                semaphore.notify();
            }
            /*
            for( int k = 0; k < fast_region.length ; k++){
                if( fast_region[k] != bit_region[k] ){
                    System.out.println(seq.Name + " mismatch: at " + k + " fa:"+ fast_region[k] + ",2bit:" + bit_region[k] );
                    TestPassed = false;
                    ++c;
                    if( c >= limit){
                        break;
                    }
                }
            }*/

        }
    }
    /**
     * the code below creates a quick and dirty pipeline to use one process to read the 
     * fa file and another process to read teh 2bit file.
     */
    
    
    static class bucket{
        public bucket( char fa[] , String fa_name){
            this.fa = fa;
            faName = fa_name;
        }
        char fa[];
        
        String faName;
    }
    static boolean done = false;
    final static Object semaphore = new String("");
    
    static LinkedList<bucket> queue = new LinkedList<bucket>();
    
    static public class worker extends Thread{
        
        @Override
        public void run(){
            System.out.println("started ");
            int count = 0;
            
            TwoBit mine;
            AnnotatedSeqGroup grp_bit = new AnnotatedSeqGroup("bgroup");
            mine = new TwoBit( twobit, "", grp_bit);
            
            while( queue.size() > 0 || !done ){
                try {
                    bucket b = null;
                    synchronized( queue){
                        b = queue.poll();
                    }
                    if( b == null){
                        if( queue.size() == 0 && done ){
                            break;
                        }
                        synchronized( semaphore){
                            try {
                                semaphore.wait();
                            } catch (InterruptedException ex) {
                                Logger.getLogger(Fast2BitCompare.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }else{
                        
                        
                        BioSeq bit_seq = map.get( b.faName );
                        SeqSpan bit_span = new SimpleSeqSpan(0 , bit_seq.getLength(), bit_seq );
                        String bit_reg;  

                        bit_reg = mine.getRegionResidues( bit_span);

                        char[] bit = bit_reg.toCharArray();

                        if( b.fa.length != bit.length ){
                            System.out.println(" Bad fasta:" + b.fa.length + " and twobit" + bit.length );
                            TestPassed = false;
                        }else{
                            //System.out.println( " The size checks out ");
                        }
                        
                        
                        for( int k = 0; k < b.fa.length ; k++){
                        if( b.fa[k] != bit[k] ){
                            System.out.println( b.faName + " mismatch: at " + k + " fa:"+ b.fa[k] + ",2bit:" + bit[k] );
                            TestPassed = false;
                        }
                        

                        if( count++ % 10000000 == 0){
                            System.out.println( this.getName() + "Sequence number " + count );
                        }
                        
                    }
                    }
            
                } catch (Exception ex) {
                        Logger.getLogger(Fast2BitCompare.class.getName()).log(Level.SEVERE, null, ex);
                }
        
            }
        }
    
    }
}