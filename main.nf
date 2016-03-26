#!/usr/bin/env nextflow

params.hml           = "${baseDir}/tutorial/ex00_ngsp_expected.xml"
params.output        = "${baseDir}/tutorial/output"
params.imgtdir       = file("${baseDir}/imgt")
params.imgt          = "3200"
imgtdb               = "${params.imgtdir}/${params.imgt}"

outputDir    = file("${params.output}")
outputDir.mkdirs()

expectedFile = file("${params.hml}")
expectedFile.copyTo("$outputDir/ex11_ngsp_expected.xml")


// Filtering out the failed subjects
process filterExpectedHml{

  input:
    set file(expected) from file("${params.hml}")

  output:
    set file(expected), file('*.gz') into fastqFiles mode flatten

  """
    ngs-extract-consensus -i ${expected}
  """
}

subjectIdFiles = fastqFiles.map{ hml, fileIn ->
  tuple(subjectId(fileIn), fileIn, hml ) 
}


//Blasting the 
process blastn{
  tag{ subject }

  input:
    set subject, file(subjectFastq), file(hmlFailed) from subjectIdFiles

  output:
    set subject, file {"${subject}.failed.txt"}  into blastObservedFile mode flatten
    set file {"${subject}.failed.txt"}  into finalFailedObserved mode flatten
    set subject, file{"${hmlFailed}"} into failedHmlFiles

  """
    gzcat ${subjectFastq} | blastn -db $imgtdb -outfmt 6 -query - > blast.out
    ngs-extract-blast -i blast.out -f ${subjectFastq} > ${subject}.failed.txt
  """
}

blastObservedSubjects = blastObservedFile
.collectFile() { subject, blast ->
       [ "${subject}.txt", blast.text ]
   }
.map{ path ->
  tuple( path, blastSubjectId(path), path) 
} 

finalFailedObserved
.collectFile() {  blast ->
       [ "ex11_ngsp_observed.txt", blast.text ]
   }
.subscribe { file -> copyToFailedDir(file) }


failedSubjects = blastObservedSubjects 
.map{ observed ->
  [observed[0], observed[1]]
}

process validateInterpretation {
  tag{ subject }

  input:
    set file(observed), subject from failedSubjects
    set file(expected) from file("${params.hml}")

  output:
    set  file("${subject}_validate.txt") into failedValidated
    stdout inputDir

  """
    ngs-extract-expected-haploids -i ${expected} | ngs-validate-interpretation -b ${observed}  > ${subject}_validate.txt
    echo $outputDir
  """

}

failedValidatedFiles = failedValidated
.collectFile() { validated ->
       [ "ex11_ngsp_validated.txt", validated.text ]
   }
failedValidatedFiles.subscribe { file -> copyToFailedDir(file) }


process generateReport {
  
  input:
    set infiles from inputDir

  """
    ngs-validation-report -t failedReport -f -v 1 -p $outputDir -d $infiles
  """ 

}


def copyToFailedDir (file) { 
  log.info "Copying ${file.name} into: $outputDir"
  file.copyTo(outputDir)
}


def subjectId(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(\d{4}-\d{4}-\d{1})_\d{1,2}_\d{1,2}.fa.gz$/
  return loc[0][1]
}

def getName(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(ex\d{1,3}).txt$/
  return name
}


def blastSubjectId(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(\d{1,4}-\d{1,4}-\d{1}).txt$/
  return loc[0][1]
}

def blastExp2(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /\d{4}-\d{4}-\d{1}.(ex\d{1,3}).txt$/
  //return loc[0][1]
  return name
}

def blastExp(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(ex\d{1,3})_ngsp_observed.txt$/
  return loc[0][1]
}









