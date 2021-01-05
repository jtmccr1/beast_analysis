#!/usr/bin/env nextflow
nextflow.enable.dsl=2                                                                                                    

params.seed=new Date().getTime()
println(params.seed)
def random= new Random(params.seed)



params.n=2
params.burnin=10;
params.thinning_factor=1;
params.beast_options="";

params.tree_burnin=params.burnin;
params.tree_thinning_factor=params.thinning_factor;


beast_seeds = [];



for(int i=0;i<params.n;i++){
    beast_seeds.add(random.nextInt() & Integer.MAX_VALUE)
}

params.outDir="./"
process beast{
    	publishDir "${params.outDir}/", mode:"copy", overwrite:"true"
        input:
               tuple path(xml_file), val(seed)
        output:
                tuple val("${xml_file.name}"), path("*log"), emit: logs
                tuple val("${xml_file.name}"), path("*trees"), emit:trees
                path("*ops")
                path("*out")
                path("*chkpt") optional true
"""
beast   -prefix ${seed}_ -seed ${seed} ${params.beast_options}  ${xml_file} > ${seed}_${xml_file.name}.out
"""
}

process combine_logs{
    publishDir "${params.outDir}/combined/", mode:"copy", overwrite:"true"
	errorStrategy 'finish'
    input:
        tuple val(xml), path(logs)
    output:
            path("${xml}.log")

"""
RESAMPLE=\$(tail -n2 ${logs[0]} | awk '{print \$1}'| sort  |paste -sd- - | bc | awk '{printf "%0.f", \$1*$params.thinning_factor*$params.n}')
BURNIN=\$(tail -n1 ${logs[0]}| awk '{printf "%.0f", \$1*($params.burnin/100)}')
logcombiner -burnin \${BURNIN} \
-resample \${RESAMPLE}  $logs  ${xml}.log
"""
}


process combine_trees{
    	publishDir "${params.outDir}/combined/", mode:"copy", overwrite:"true"

	 errorStrategy 'finish'
        input:
        tuple val(xml), path(trees)
        output:
            path("${xml}.trees")
"""
RESAMPLE=\$(awk '/^tree/{split(\$2,a,"_"); printf "%s\\n",a[2]}' ${trees[0]} | tail -n2 | sort  |paste -sd- - | bc| awk '{printf "%0.f", \$1*$params.tree_thinning_factor*$params.n}')
BURNIN=\$(awk '/^tree/{split(\$2,a,"_"); printf "%s\\n",a[2]}' ${trees[0]} | tail -n1| awk '{printf "%.0f", \$1*($params.tree_burnin/100)}')
logcombiner  -trees   -burnin \${BURNIN} \
-resample \${RESAMPLE}  $trees  ${xml}.trees
"""

}
process mcc{
    	publishDir "${params.outDir}/combined/", mode:"copy", overwrite:"true"
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
