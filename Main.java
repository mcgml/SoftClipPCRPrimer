package nhs.genetics.cardiff;

import htsjdk.samtools.*;
import nhs.genetics.cardiff.framework.BEDFile;
import nhs.genetics.cardiff.framework.Clip;
import nhs.genetics.cardiff.framework.GenomicLocation;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Programme for soft clipping PCR primers based on start-stop locations
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2016-10-13
 */
public class Main {
    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final String program = "SoftClipPCRPrimers";
    private static final String version = "1.0.0";

    public static void main(String[] args) throws IOException {

        CommandLineParser commandLineParser = new DefaultParser();
        CommandLine commandLine = null;
        HelpFormatter formatter = new HelpFormatter();
        Options options = new Options();

        options.addOption("I", "Input", true, "Path to input BAM file");
        options.addOption("T", "Targets", true, "Path to BED file");
        options.addOption("O", "Output", true, "Path to output BAM file");

        try {
            commandLine = commandLineParser.parse(options, args);

            if (!commandLine.hasOption("I") || !commandLine.hasOption("T") || ! commandLine.hasOption("O")){
                throw new NullPointerException("Incorrect arguments");
            }

        } catch (ParseException | NullPointerException e){
            formatter.printHelp(program + " " + version, options);
            log.log(Level.SEVERE, e.getMessage());
            System.exit(-1);
        }

        File inputSamOrBamFile = new File(commandLine.getOptionValue("I"));
        File outputSamOrBamFile = new File(commandLine.getOptionValue("O"));
        File bedFile = new File(commandLine.getOptionValue("T"));
        ArrayList<GenomicLocation> genomicLocations = null;

        log.log(Level.INFO, "Reading BED file: " + bedFile + " ...");
        try {
            genomicLocations = BEDFile.getBedFeatures(bedFile);
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not read BED file: " + e.getMessage());
            System.exit(-1);
        }

        log.log(Level.INFO, "Reading BAM file: " + inputSamOrBamFile + " ...");
        try (SamReader samReader = SamReaderFactory.makeDefault().open(inputSamOrBamFile)){

            SAMFileHeader samFileHeader = samReader.getFileHeader();
            samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            SAMProgramRecord samProgramRecord = samFileHeader.createProgramRecord();
            samProgramRecord.setCommandLine(String.join(" ", args));
            samProgramRecord.setProgramName(program);
            samProgramRecord.setProgramVersion(version);

            log.log(Level.INFO, "Processing reads, writing to " + outputSamOrBamFile.getName() + " ...");
            try (SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, outputSamOrBamFile)){

                //loop over BED records
                for (GenomicLocation genomicLocation : genomicLocations){
                    log.log(Level.INFO, "Inspecting region: " + genomicLocation + " ...");

                    //query alignments
                    SAMRecordIterator samRecordIterator = samReader.queryAlignmentStart(genomicLocation.getContig(), genomicLocation.getStartPosition());

                    samRecordIterator.stream()
                            .filter(samRecord -> !samRecord.getReadUnmappedFlag())
                            .filter(samRecord -> !samRecord.getNotPrimaryAlignmentFlag())
                            .filter(samRecord -> !samRecord.getSupplementaryAlignmentFlag())
                            .filter(samRecord -> !samRecord.getReadFailsVendorQualityCheckFlag())
                            .filter(samRecord -> samRecord.getAlignmentEnd() == genomicLocation.getEndPosition())
                            .filter(samRecord -> {
                                return samRecord.getCigar().getFirstCigarElement().getLength() >= genomicLocation.getUpstreamPrimerLength() &&
                                        samRecord.getCigar().getFirstCigarElement().getOperator().equals(CigarOperator.M) &&
                                        samRecord.getCigar().getLastCigarElement().getLength() >= genomicLocation.getDownstreamPrimerLength() &&
                                        samRecord.getCigar().getLastCigarElement().getOperator().equals(CigarOperator.M);
                            })
                            .forEach(samRecord -> {
                                samFileWriter.addAlignment(Clip.softClip(samRecord, genomicLocation.getUpstreamPrimerLength(), genomicLocation.getDownstreamPrimerLength()));
                            });

                    samRecordIterator.close();

                }

            }

        }

    }

}