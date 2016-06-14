package nhs.genetics.cardiff;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import static java.lang.Math.toIntExact;

/**
 * Created by ml on 19/04/2016.
 */
public class BedRecord {

    private final static BEDCodec bedCodec = new BEDCodec(BEDCodec.StartOffset.ONE);
    private BEDFeature bedFeature;
    private int downstreamPrimerLength, upstreamPrimerLength, thickStart, thickEnd;

    public BedRecord(String line){
        bedFeature = bedCodec.decode(line);

        String[] fields = line.split("\t");

        thickStart = Integer.parseInt(fields[6]) + 1; //1-based
        thickEnd =Integer.parseInt(fields[7]);
    }

    public BEDFeature getBedFeature() {
        return bedFeature;
    }
    public int getThickStart() {
        return thickStart;
    }
    public int getThickEnd() {
        return thickEnd;
    }
    public int getDownstreamPrimerLength() {
        return bedFeature.getEnd() - thickEnd;
    }
    public int getUpstreamPrimerLength() {
        return thickStart - bedFeature.getStart();
    }
    public void setThickStart(int thickStart) {
        this.thickStart = thickStart;
    }
    public void setThickEnd(int thickEnd) {
        this.thickEnd = thickEnd;
    }

}
