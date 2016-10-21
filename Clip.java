package nhs.genetics.cardiff.framework;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CigarUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Class for parsing clipping cigars
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2016-10-14
 */
public class Clip {
    public static SAMRecord softClip(SAMRecord samRecord, int clipFromStart, int clipFromEnd){
        int clipFrom;

        if (clipFromStart > 0) {
            clipFrom = (samRecord.getReadLength() - clipFromStart) + 1;

            samRecord.setAlignmentStart(samRecord.getUnclippedStart() + clipFromStart);
            samRecord.setCigar(new Cigar(invertCigarElements(CigarUtil.softClipEndOfRead(clipFrom, invertCigarElements(samRecord.getCigar().getCigarElements())))));
        }

        if (clipFromEnd > 0){
            clipFrom = (samRecord.getReadLength() - clipFromEnd) + 1;

            samRecord.setCigar(new Cigar(CigarUtil.softClipEndOfRead(clipFrom, samRecord.getCigar().getCigarElements())));
        }
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
