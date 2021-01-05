#!/usr/bin/env nextflow
nextflow.enable.dsl=2                                                                                                    

params.seed=new Date().getTime()
println(params.seed)
def random= new Random(params.seed)



params.n=2
params.burnin=[10,20];
params.thinning_factor=[2,4];
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
        tuple val(xml), path(logs), val(burnin), val(thinning)
    output:
            path("${xml}.b${burnin}.thin${thinning}.log")

"""
RESAMPLE=\$(tail -n2 ${logs[0]} | awk '{print \$1}'| sort  |paste -sd- - | bc | awk '{printf "%0.f", \$1*$thinning}')
BURNIN=\$(tail -n1 ${logs[0]}| awk '{printf "%.0f", \$1*($burnin/100)}')
logcombiner -burnin \${BURNIN} \
-resample \${RESAMPLE}  $logs  ${xml}.b${burnin}.thin${thinning}.log
"""
}


process combine_trees{
    	publishDir "${params.outDir}/combined/", mode:"copy", overwrite:"true"

	 errorStrategy 'finish'
        input:
        tuple val(xml), path(trees),val(burnin), val(thinning)
        output:
            path("${xml}.b${burnin}.thin${thinning}.trees")
"""
RESAMPLE=\$(awk '/^tree/{split(\$2,a,"_"); printf "%s\\n",a[2]}' ${trees[0]} | tail -n2 | sort  |paste -sd- - | bc| awk '{printf "%0.f", \$1*$thinning}')
BURNIN=\$(awk '/^tree/{split(\$2,a,"_"); printf "%s\\n",a[2]}' ${trees[0]} | tail -n1| awk '{printf "%.0f", \$1*($burnin/100)}')
logcombiner  -trees   -burnin \${BURNIN} \
-resample \${RESAMPLE}  $trees  ${xml}.b${burnin}.thin${thinning}.trees
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
burnin_ch=channel.from(params.burnin)
thinning_ch=channel.from(params.thinning_factor)
tree_burnin_ch=channel.from(params.tree_burnin)
tree_thinning_ch=channel.from(params.tree_thinning_factor)

workflow {
    
    beast(xml_ch.combine(channel.from(beast_seeds)))
    combine_logs(beast.out.logs.groupTuple(size:params.n).combine(burnin_ch).combine(thinning_ch))
   	combine_trees(beast.out.trees.groupTuple(size:params.n).combine(tree_burnin_ch).combine(tree_thinning_ch))
	mcc(combine_trees.out)
}
