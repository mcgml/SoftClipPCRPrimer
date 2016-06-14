package nhs.genetics.cardiff;

import htsjdk.samtools.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {
    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) throws IOException {

        if (args.length != 2) {
            System.err.println("Usage: <bam> <bed>");
            System.err.println("Bam file must be name sorted and processed using fixmates.");
            System.err.println("Bed file must contain amplicons. Target regions must be specified using thick and thin regions");
            System.exit(1);
        }

        int numberOfReads = 0;
        File inputSamOrBamFile = new File(args[0]);
        File outputSamOrBamFile = new File(args[0].replace(".bam", "_clipped.bam"));
        SAMRecord firstInPairRecord = null, secondInPairRecord = null;
        HashMap<String, Integer> ampliconCounts = new HashMap<>();

        //prepare bam reader
        SamReader samReader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
        SAMFileHeader samFileHeader = samReader.getFileHeader();

        if (!samFileHeader.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
            throw new IllegalArgumentException("Supplied bam must be name sorted");
        }

        //prepare writer
        samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, outputSamOrBamFile);

        //read BED file and store in mem
        HashMap<String, BedRecord> bedRecords = getBedFeatures(new File(args[1]));

        //init count map
        for (Map.Entry<String, BedRecord> amplicon : bedRecords.entrySet()){
            ampliconCounts.put(amplicon.getValue().getBedFeature().getName(), 0);
        }

        //loop over sam/bam
        log.log(Level.INFO, "Processing reads...");
        for (SAMRecord samRecord : samReader){
            numberOfReads++;

            //skip unmapped reads and read pairs mapping to the same strand
            if (samRecord.getReadUnmappedFlag() || samRecord.getMateUnmappedFlag() || samRecord.getNotPrimaryAlignmentFlag()){
                continue;
            }

            if (samRecord.getFirstOfPairFlag()){
                firstInPairRecord = samRecord;
            } else if (samRecord.getSecondOfPairFlag()) {
                secondInPairRecord = samRecord;

                //check pairing - sanity check
                if (!firstInPairRecord.getReadName().equals(secondInPairRecord.getReadName())) {
                    throw new IllegalArgumentException("Unexpected pairing of reads: " + firstInPairRecord.getReadName() + " " + secondInPairRecord.getReadName());
                }

                //soft clip primers
                if (firstInPairRecord.getProperPairFlag() && secondInPairRecord.getProperPairFlag()){
                   String coor = firstInPairRecord.getContig() + ":" + firstInPairRecord.getUnclippedStart() + "-" + secondInPairRecord.getUnclippedEnd();

                    if (bedRecords.containsKey(coor)){

                        //skip reads with ins/del in primer
                        if (firstInPairRecord.getCigar().getFirstCigarElement().getOperator() != CigarOperator.M) continue;
                        if (secondInPairRecord.getCigar().getFirstCigarElement().getOperator() != CigarOperator.M) continue;
                        if (firstInPairRecord.getCigar().getFirstCigarElement().getLength() < bedRecords.get(coor).getUpstreamPrimerLength()) continue;
                        if (secondInPairRecord.getCigar().getFirstCigarElement().getLength() < bedRecords.get(coor).getDownstreamPrimerLength()) continue;

                        //count amplicon frequencies
                        ampliconCounts.put(bedRecords.get(coor).getBedFeature().getName(), ampliconCounts.get(bedRecords.get(coor).getBedFeature().getName()) + 1);

                        //R1 upstream & downstream oligos
                        firstInPairRecord = SoftClipRead.softClipFivePrime(firstInPairRecord, bedRecords.get(coor).getUpstreamPrimerLength());
                        firstInPairRecord = SoftClipRead.softClipThreePrime(firstInPairRecord, firstInPairRecord.getUnclippedEnd() - bedRecords.get(coor).getThickEnd());

                        //R2 downstream & upstream oligos
                        secondInPairRecord = SoftClipRead.softClipThreePrime(secondInPairRecord, bedRecords.get(coor).getDownstreamPrimerLength());
                        secondInPairRecord = SoftClipRead.softClipFivePrime(secondInPairRecord, bedRecords.get(coor).getThickStart() - secondInPairRecord.getUnclippedStart());

                        samFileWriter.addAlignment(firstInPairRecord);
                        samFileWriter.addAlignment(secondInPairRecord);
                    }

                }

            } else {
                throw new IllegalArgumentException("Read is neither first or second in pair: " + samRecord.getReadName());
            }

        }

        log.log(Level.INFO, "Processed " + numberOfReads + " reads");

        //print per amplicon frequencies
        for (Map.Entry<String, Integer> amplicon : ampliconCounts.entrySet()){
            System.out.println(amplicon.getKey() + "\t" + amplicon.getValue());
        }

        samReader.close();
        samFileWriter.close();
    }

    private static HashMap<String, BedRecord> getBedFeatures(File filePath) throws IOException {
        HashMap<String, BedRecord> bedRecords = new HashMap<>();
        String line;

        try (BufferedReader bedReader = new BufferedReader(new FileReader(filePath))){
            while((line = bedReader.readLine())!= null) {

                if (!line.equals("")){

                    BedRecord bedRecord = new BedRecord(line);
                    String coor = bedRecord.getBedFeature().getContig() + ":" + bedRecord.getBedFeature().getStart() + "-" + bedRecord.getBedFeature().getEnd();

                    if (bedRecords.containsKey(coor)){
                        BedRecord previousBedRecord = bedRecords.get(coor);

                        //shrink thick region to smallest size
                        if (bedRecord.getThickStart() > previousBedRecord.getThickStart()){
                            previousBedRecord.setThickStart(bedRecord.getThickStart());
                        }
                        if (bedRecord.getThickEnd() < previousBedRecord.getThickEnd()){
                            previousBedRecord.setThickEnd(bedRecord.getThickEnd());
                        }

                        bedRecords.put(coor, previousBedRecord);

                    } else {
                        bedRecords.put(coor, bedRecord);
                    }

                }

            }

            bedReader.close();
        }

        return  bedRecords;
    }


}