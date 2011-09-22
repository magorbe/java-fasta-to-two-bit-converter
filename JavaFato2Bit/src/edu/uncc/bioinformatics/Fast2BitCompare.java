/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.uncc.bioinformatics;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.SeqSymmetry;
import com.affymetrix.genometryImpl.span.SimpleSeqSpan;
import com.affymetrix.genometryImpl.symloader. Fasta;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import java.io.File;
import java.net.URI;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jfvillal
 */
public class Fast2BitCompare {
	public static boolean TestPassed = true;
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            String[] files = args[0].split("@");
            HashMap<String,String> map = new HashMap<String,String>();
            for( int i = 0; i < files.length; i++){
                File f = new File( files[i]);
                map.put(f.getName().replace(".fasta", ""), files[i]);    
            }
            
            URI twobit = URI.create( args[1]);
            AnnotatedSeqGroup grp_bit = new AnnotatedSeqGroup("bgroup");
            TwoBit bit = new TwoBit( twobit, "", grp_bit);
            
            List<BioSeq> bit_lst = bit.getChromosomeList();
            for( BioSeq seq : bit_lst ){
                
                String name = map.get(  seq.getID() );
                if( name == null){
                    System.out.println("sequence not found in fasta. Validatation failed");
                    TestPassed = false;
                    break;
                }
                URI res  = URI.create( name );
                AnnotatedSeqGroup grp_fas = new AnnotatedSeqGroup("agroup");
                Fasta fast = new Fasta( res, "test" , grp_fas  );                

                List<BioSeq> lst = fast.getChromosomeList();
                SeqSpan span = new SimpleSeqSpan(0, lst.get(0).getLength() , lst.get(0));
                String fast_region = fast.getRegionResidues( span );
                
                
                
                
                SeqSpan bit_span = new SimpleSeqSpan(0 , seq.getLength(), seq );
                String bit_region = bit.getRegionResidues( bit_span);

                if( fast_region.length() != bit_region.length()){
                    System.out.println(" Bad fasta:" + fast_region.length() + " and twobit" + bit_region.length() );
                    TestPassed = false;
                }else{
                    System.out.println( " The size checks out ");
                }

                int limit = 100;
                int c = 0;
                for( int i = 0; i < fast_region.length() ; i++){
                    if( fast_region.charAt(i) != bit_region.charAt(i) ){
                        System.out.println("One does not checkout " + fast_region.charAt(i) + "," + bit_region.charAt(i) );
                        TestPassed = false;
                        ++c;
                        if( c >= limit){
                            break;
                        }
                    }
                }
            }
            System.out.println("Done");
            
        } catch (Exception ex) {
            Logger.getLogger(Fast2BitCompare.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
}
