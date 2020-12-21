#!/usr/bin/env nextflow
nextflow.enable.dsl=2                                                                                                    

params.seed=new Date().getTime()
println(params.seed)
def random= new Random(params.seed)



params.n=1
params.chainLength=1000000  
params.burnin=10
params.xml=""
params.beast_options=""
beast_seeds = []



for(int i=0;i<params.n;i++){
    beast_seeds.add(random.nextInt() & Integer.MAX_VALUE)
}

params.outDir="./"
process beast{
    	publishDir "${params.outDir}/${workflow.runName}/", mode:"copy", overwrite:"true", pattern: "*log"
    	publishDir "${params.outDir}/${workflow.runName}/", mode:"copy", overwrite:"true", pattern: "*out"
    	publishDir "${params.outDir}/${workflow.runName}/", mode:"copy", overwrite:"true", pattern: "*trees"
    	publishDir "${params.outDir}/${workflow.runName}/", mode:"copy", overwrite:"true", pattern:"*ops"
        input:
               tuple path(xml_file), val(seed)
        output:
                tuple val("${xml_file.name}"), path("*log"), emit: logs
                tuple val("${xml_file.name}"), path("*trees"), emit:trees
                path("*ops")
                path("*out")
"""
beast   -prefix ${seed}_ -seed ${seed} ${params.beast_options}  ${xml_file} > ${seed}_${xml_file.name}.out
"""


}


process combine_logs{
    	publishDir "${params.outDir}/${workflow.runName}/combined/", mode:"copy", overwrite:"true"
	 errorStrategy 'finish'
        input:
        tuple val(xml), path(logs)
        output:
            path("${xml}.log")

"""
logcombiner   -burnin ${(params.chainLength/params.burnin).round(0)} \
-resample ${(params.chainLength/10000)*params.n}  $logs  ${xml}.log
"""
}
process combine_trees{
    	publishDir "${params.outDir}/${workflow.runName}/combined/", mode:"copy", overwrite:"true"

	 errorStrategy 'finish'
        input:
        input:
        tuple val(xml), path(trees)
        output:
            path("${xml}.trees")

"""
logcombiner  -trees -burnin ${(params.chainLength/params.burnin).round(0)} \
-resample ${(params.chainLength/10000)*params.n}  $trees  ${xml}.trees
"""
}
process mcc{
    	publishDir "${params.outDir}/${workflow.runName}/combined/", mode:"copy", overwrite:"true"
	 errorStrategy 'finish'
        input:
	input:
	    path(tree)
	output:
	    path("${tree.name}.mcc.tree")
"""
treeannotator  $tree ${tree.name}.mcc.tree
"""
}


xml_ch=channel.fromPath(params.xml).flatMap()

workflow {
    
    beast(xml_ch.combine(channel.from(beast_seeds)))
    combine_logs( beast.out.logs.groupTuple(size:params.n))
    combine_trees(beast.out.trees.groupTuple(size:params.n))
    mcc(combine_trees.out)
}
