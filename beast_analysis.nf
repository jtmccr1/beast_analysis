#!/usr/bin/env nextflow
nextflow.enable.dsl=2                                                                                                    

def random= new Random(params.seed)

beast_seeds=[];
def outputDateFormat = "yyyy-MM-dd'T'HH:mm:ssZ"
now = new Date().format(outputDateFormat)
for(int i=0;i<params.n;i++){
    beast_seeds.add(random.nextInt() & Integer.MAX_VALUE)
}

params.outDir="./"

process make_config {
publishDir "${params.outDir}/", mode:"copy", overwrite:"true"
output:
	path(("*.params.txt"))
"""
echo "params {
    n=${params.n};
    burnin=${params.burnin};
    thinning_factor=${params.thinning_factor};
    beast_options=\\"${params.beast_options}\\";
    jvm_options=\\"${params.jvm_options}\\";
    seed=${params.seed};
    save_every=${params.save_every};
    tree_burnin=${params.tree_burnin};
    tree_thinning_factor=${params.tree_thinning_factor};
    xml = ${params.xml};
    data = ${params.data};
    template=${params.template};
    process_beast=${params.process_beast};
}
" > ${now}.params.txt

"""

}

process beastgen{
	    publishDir "${params.outDir}/", mode:"copy", overwrite:"true"
		input:
			tuple path(data), path(template)
		output:
			path("*xml")
			
"""
cp $template ./local_template
beastgen -date_order -1 -date_prefix "|" -date_precision -D "outputStem=${data.name.take(data.name.lastIndexOf('.'))}.${template.name.take(template.name.lastIndexOf('.'))}"	local_template $data ${data.name.take(data.name.lastIndexOf('.'))}.${template.name.take(template.name.lastIndexOf('.'))}.xml
"""
}
process beast{
    	publishDir "${params.outDir}/", mode:"copy", overwrite:"true"
        input:
               tuple path(xml_file), val(seed)
        output:
                tuple val("${xml_file.name.take(xml_file.name.lastIndexOf('.'))}"), path("*log"), emit: logs
                tuple val("${xml_file.name.take(xml_file.name.lastIndexOf('.'))}"), path("*trees"), emit:trees
                path("*ops")
                path("*out")
                path("*chkpt") optional true
"""
export JAVA_TOOL_OPTIONS='${params.jvm_options}'
beast  ${(params.save_every>0? "-save_every ${params.save_every} -save_state ${seed}_${xml_file.name.take(xml_file.name.lastIndexOf('.'))}.chkpt":'')}  -prefix ${seed}_ -seed ${seed} ${params.beast_options}  ${xml_file} > ${seed}_${xml_file.name.take(xml_file.name.lastIndexOf('.'))}.out
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


burnin_ch=channel.from(params.burnin)
thinning_ch=channel.from(params.thinning_factor)
tree_burnin_ch=channel.from(params.tree_burnin)
tree_thinning_ch=channel.from(params.tree_thinning_factor)

if(params.xml==null &&(params.data==null || params.template==null)){
			throw new Exception("must provide either xml or datafile and template")
		}

workflow {
        make_config()
        println("running beast workflow with seed ${params.seed}")
        if(params.xml!=null){
            xml_ch=channel.fromPath(params.xml).flatMap()
            beast(xml_ch.combine(channel.from(beast_seeds)))
        }else{
            data_ch=channel.fromPath(params.data).flatMap()
            template_ch=channel.fromPath(params.template).flatMap()
            beastgen(data_ch.combine(template_ch))
            beast(beastgen.out.combine(channel.from(beast_seeds)))
        }
    if(params.process_beast){
        if(params.n>1){
            combine_logs(beast.out.logs.groupTuple(size:params.n).combine(burnin_ch).combine(thinning_ch))
            combine_trees(beast.out.trees.groupTuple(size:params.n).combine(tree_burnin_ch).combine(tree_thinning_ch))
            mcc(combine_trees.out)
        }
    }
}



