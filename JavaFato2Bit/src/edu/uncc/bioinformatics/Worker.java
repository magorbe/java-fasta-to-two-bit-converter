package edu.uncc.bioinformatics;

import java.net.URI;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.span.SimpleSeqSpan;
import com.affymetrix.genometryImpl.symloader.TwoBit;

public class Worker extends Thread{
	URI Url;
	ArrayBlockingQueue<Bucket> Queue ;
    public Worker( URI url ,     ArrayBlockingQueue<Bucket> queue  ){
    	Url = url;
    	Queue = queue;
    }
    public boolean getPassed(){
    	return TestPassed;
    }
    static boolean TestPassed;
    @Override
    public void run(){
    	try {
            long count = 0;
            
            TwoBit mine;
            AnnotatedSeqGroup grp_bit = new AnnotatedSeqGroup("bgroup");
            mine = new TwoBit( Url, "", grp_bit);
            
            HashMap<String,BioSeq> map = new HashMap<String,BioSeq>();
            
            List<BioSeq> bit_lst = mine.getChromosomeList();
			
            for( BioSeq seq : bit_lst ){
                map.put(seq.getID(), seq);    
            }
            bit_lst = null;
            int num = 1;
            int total = map.size();
            while(true ){
                try {
                    Bucket b = null;
                    
                    b = Queue.poll();
                    
                    if( b == null){
                        
                    	Thread.sleep(100);
                        
                    }else{
                    	
                        if( b.new_seq ){
                        	//System.out.println( num++ + " of " + total );
                        	if( b.getDone() ){
                        		System.out.println("Done");
                        		break;
                        	}
                        	continue;
                        }
                       /*if( b.new_seq ){
                    	   	grp_bit = new AnnotatedSeqGroup("bgroup");
           	            	mine = new TwoBit( twobit, "", grp_bit);
           	            	continue;
                       }*/
                        BioSeq bit_seq = map.get( b.faName );
                        SeqSpan bit_span = new SimpleSeqSpan(0 , bit_seq.getLength(), bit_seq );
                        String bit_reg;  

                        bit_reg = mine.getRegionResidues( bit_span );

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
                        

                        if( count++ % 100000000 == 0){
                            System.out.print("."  );
                        }
                       
                        }
                        
                        //free up memory
                        map.remove(b.faName);
                    }
            
                } catch (Exception ex) {
                        Logger.getLogger(Fast2BitCompare.class.getName()).log(Level.SEVERE, null, ex);
                }
        
            }
            System.out.println();
    	} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

}
