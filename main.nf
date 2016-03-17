#!/usr/bin/env nextflow

params.expected      = "${baseDir}/ex00_ngsp_expected.xml"
//params.finalDir    = "${baseDir}/tutorial"
params.failed        = "${baseDir}/tutorial/failed"
params.experiment    = "ex00"
//params.resolution  = 2
params.imgtdb        = "${baseDir}/imgt/3200"

params.validatedDir = "${baseDir}/validated"
params.testExpected = "${baseDir}/ex00_ngsp_expected.xml,${baseDir}/ex01_ngsp_expected.xml"
expectedFiles       = Channel.from( "${params.testExpected}" )


validatedDir = file("${params.validatedDir}")
failedDir    = file("${params.failed}")

failedDir.mkdirs()
validatedDir.mkdirs()

failedSubjects = Channel.fromPath('failed_subjects.csv').splitCsv()
.collectFile() { row ->
       [ "${row[0]}.txt", row[1] + '\n' ]
   }
.map{path -> 
  tuple(exp(path), path )
}blastObservedFile


// Filtering out the failed subjects
process filterExpectedHml{
  tag{ exp }

  input:
    set exp, file(failedFile) from failedSubjects

  output:
    set exp, file{"${exp}_ngsp_expected.xml"}, file('*.gz') into fastqFiles mode flatten
    file{"${exp}_ngsp_expected.xml"} into failedExpectedFiles

  """
    ngs-filter-samples -i ${params.expected} -s ${failedFile} > ${exp}_ngsp_expected.xml
    ngs-extract-consensus -i ${exp}_ngsp_expected.xml
  """
}

subjectIdFiles = fastqFiles.map{ exp, hml, fileIn ->
  tuple(exp, subjectId(fileIn), fileIn, hml ) 
}

//Copying the new expected files
failedExpectedFiles
.subscribe { file -> copyToFailedDir(file) }


//Blasting the 
process blastn{
  tag{ subject }

  input:
    set exp, subject, file(subjectFastq), file(hmlFailed) from subjectIdFiles

  output:
    set exp, subject, file {"${subject}.failed.txt"}  into blastObservedFile mode flatten
    set exp, file {"${subject}.failed.txt"}  into finalFailedObserved mode flatten
    set exp, subject, file{"${hmlFailed}"} into failedHmlFiles

  """
    gzcat ${subjectFastq} | blastn -db ${params.imgtdb} -outfmt 6 -query - > blast.out
    ngs-extract-blast -i blast.out -f ${subjectFastq} > ${subject}.failed.txt
  """
}




blastObservedSubjects = blastObservedFile
.collectFile() { exp, subject, blast ->
       [ "${subject}.${exp}.txt", blast.text ]
   }
.map{ path ->
  tuple( blastExp2(path), blastSubjectId(path), path) 
}

finalFailedObserved
.collectFile() { exp, blast ->
       [ "${exp}_ngsp_observed.txt", blast.text ]
   }
.subscribe { file -> copyToFailedDir(file) }


failedSubjects = blastObservedSubjects 
.phase(failedHmlFiles)
.map{ observed, hml ->
  [observed[0], observed[1], observed[2], hml[2]]
}


process validateInterpretation {
  tag{ subject  }

  input:
    set exp, subject, file(observed), file(expected) from failedSubjects

  output:
    set exp, file {"${exp}_${subject}_validate.txt"} into failedValidated

  """
    ngs-extract-expected-haploids -i ${expected} | ngs-validate-interpretation -b ${observed}  > "${exp}_${subject}_validate.txt"
  """

}


failedValidatedFiles = failedValidated
.collectFile() { exp, validated ->
       [ "${exp}_ngsp_validated.txt", validated.text ]
   }
failedValidatedFiles.subscribe { file -> copyToFailedDir(file) }


failedReportCmdline = """  ngs-validation-report -d $failedDir -t failedReport -f -v 1 """
result = failedReportCmdline.execute().text



def copyToFailedDir (file) { 
  log.info "Copying ${file.name} into: $failedDir"
  file.copyTo(failedDir)
}


def subjectId(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(\d{4}-\d{4}-\d{1})_\d{1,2}_\d{1,2}.fa.gz$/
  return loc[0][1]
}

def exp(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(ex\d{1,3}).txt$/
  return loc[0][1]
}


def blastSubjectId(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(\d{1,4}-\d{1,4}-\d{1}).ex\d{1,3}.txt$/
  return loc[0][1]
}

def blastExp2(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /\d{4}-\d{4}-\d{1}.(ex\d{1,3}).txt$/
  return loc[0][1]
}

def blastExp(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /(ex\d{1,3})_ngsp_observed.txt$/
  return loc[0][1]
}









