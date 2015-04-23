package org.broadinstitute.sting.queue.qscripts.examples
import org.broadinstitute.sting.gatk.downsampling.DownsampleType
import org.broadinstitute.sting.queue.QScript
import scala.sys.process._


import org.broadinstitute.sting.queue.extensions.gatk._

class HaplotypeQueue extends QScript {

qscript =>

// Required arguments.  All initialized to empty values.
@Input(doc="The reference file for the bam files.", shortName="R")
var referenceFile: File = _ 

@Input(doc="Bam file to genotype.", shortName="I")
var bamFiles: List[File] = Nil

@Input(doc="An optional file with a list of intervals to proccess.", shortName="LF", required=false)
   var intervals: File = _

@Argument(doc="An optional with a list of command-line intervals.", shortName="L", required=false)
   var intervalsString: List[String] = Nil

@Argument(doc="An optional value to to downsample to.", shortName="dcov", required=false)
   var downsample_to_coverage : Int = _

@Argument(doc="An optional value of downsampling type.", shortName="dt", required=false)
   var downsampling_type : String = "NONE"

   def getDT : DownsampleType = {
     if (downsampling_type == "NONE")
       DownsampleType.NONE
   else if (downsampling_type == "BY_SAMPLE")
    DownsampleType.BY_SAMPLE
   else
     DownsampleType.ALL_READS
 }

@Argument(doc="An optional value to to downsample to.", shortName="stand_call_conf", required=false)
   var standard_min_confidence_threshold_for_calling : Double = _

@Argument(doc="An optional value to to downsample to.", shortName="stand_emit_conf", required=false)
   var standard_min_confidence_threshold_for_emitting : Double = _

@Input(doc="An optional file for dbsnp", shortName="D", required=false)
   var dbsnp: File = _

@Argument(doc="An optional boolean for alternate allele annotation", shortName="rdh", required=false)
   var recoverDanglingHeads: Boolean = false 

@Argument(doc="An optional boolean for alternate allele annotation", shortName="dontUseSoftClippedBases", required=false)
   var dontUseSoftClippedBases: Boolean = false 

@Argument(doc="One or more specific annotations to apply to variant calls", required=false)
   var annotation: List[String] = List("ClippingRankSumTest", "DepthPerSampleHC") 

@Argument(doc="Read filter", required=false)
  var read_filter: List[String] = Nil

@Argument(doc="scatter parameter", shortName="P")
var scatter:  Int = _


@Argument(doc="nct parameter", shortName="C")
var nct:  Int = _

@Argument(doc="Output file name", shortName="O")
var outFile: File = _

trait HCArguments extends CommandLineGATK {
        this.reference_sequence = qscript.referenceFile
        this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
        this.intervalsString = qscript.intervalsString
        this.downsample_to_coverage = qscript.downsample_to_coverage
        this.downsampling_type = qscript.getDT
        this.read_filter = qscript.read_filter
        //this.memoryLimit = 6 
}

def script() {
        val genotyper = new HaplotypeCaller with HCArguments
        genotyper.scatterCount = qscript.scatter
        genotyper.dbsnp = qscript.dbsnp
        genotyper.annotation = qscript.annotation
        genotyper.out = qscript.outFile 
        genotyper.nct = qscript.nct
	genotyper.recoverDanglingHeads = qscript.recoverDanglingHeads
	genotyper.dontUseSoftClippedBases = qscript.dontUseSoftClippedBases
        genotyper.standard_min_confidence_threshold_for_calling = qscript.standard_min_confidence_threshold_for_calling
        genotyper.standard_min_confidence_threshold_for_emitting = qscript.standard_min_confidence_threshold_for_emitting  
      for(bamFile <- bamFiles){
          genotyper.input_file :+= bamFile
        }
        add(genotyper)
}
}

