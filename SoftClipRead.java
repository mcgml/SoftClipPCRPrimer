package nhs.genetics.cardiff;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CigarUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by ml on 27/04/2016.
 */
public class SoftClipRead {

    //NB assumes no ins/del with referenceLength

    public static SAMRecord softClipFivePrime(SAMRecord samRecord, int referenceLength){
        if (referenceLength < 1) return samRecord;

        int clipFrom = (samRecord.getReadLength() - referenceLength) + 1;

        samRecord.setAlignmentStart(samRecord.getUnclippedStart() + referenceLength);
        samRecord.setCigar(new Cigar(invertCigarElements(CigarUtil.softClipEndOfRead(clipFrom, invertCigarElements(samRecord.getCigar().getCigarElements())))));

        return samRecord;
    }

    public static SAMRecord softClipThreePrime(SAMRecord samRecord, int referenceLength){
        if (referenceLength < 1) return samRecord;

        int clipFrom = (samRecord.getReadLength() - referenceLength) + 1;

        samRecord.setCigar(new Cigar(CigarUtil.softClipEndOfRead(clipFrom, samRecord.getCigar().getCigarElements())));

        return samRecord;
    }

    private static List<CigarElement> invertCigarElements(final List<CigarElement> cigar) {
        List<CigarElement> retVal = new ArrayList<>();

        for (int i = cigar.size() - 1; i >= 0; i--) {
            retVal.add(cigar.get(i));
        }

        return retVal;
    }

}
