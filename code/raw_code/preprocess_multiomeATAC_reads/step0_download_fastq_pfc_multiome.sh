## the data transfer through BaDoi's google drive

DATADIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/ATAC

## ATAC PFC data
mkdir -p ${DATADIR}; cd ${DATADIR}

# Nairobi 1
gdrive download --recursive 1NbhDwIPAPhMHYbeY3ol7IV4LkPCJKGyB
# Nairobi 2
gdrive download --recursive 1A6LrXJj0rOkX9EkylRr4l8WrFvlV49Tm
# Memphis 2
gdrive download --recursive 1Co-tYSTnBRkFdA-xV8zr2rG6SeY-ZQkF







# ## ATAC for monkey London
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/ATAC/London
# LINKFILE='https://drive.google.com/file/d/16kLYz56gwEsO9QFtEC2lSPT5lsXNJu7y/view?usp=sharing, https://drive.google.com/file/d/1De197TTWH6vEGhj2dGd-x47CIBZGftIM/view?usp=sharing, https://drive.google.com/file/d/1Fzquxm9ixMr_paZ1qs_87ZSfIO2myvTv/view?usp=sharing, https://drive.google.com/file/d/1JVsI-2yHNbgKKjyzz6zY27XLcOun_hh-/view?usp=sharing, https://drive.google.com/file/d/1KDUS-KJuzj47HaH_iYSh9oJkTqGrGSA6/view?usp=sharing, https://drive.google.com/file/d/1L3jc65aSnME-bnLnForhyiLoeFW7dSpi/view?usp=sharing, https://drive.google.com/file/d/1NdQ9BaDeg6OtOvlHuWD8sab3wbC5N4Gf/view?usp=sharing, https://drive.google.com/file/d/1RR9V1K9p41XAFOMjLCNzxeN1mNEngOgo/view?usp=sharing, https://drive.google.com/file/d/1WKzjQTSyb0izhCwKKOYQTwJzQavKxGQ0/view?usp=sharing, https://drive.google.com/file/d/1YiiSEVoS8i_STRb7HrWG9kOa0Qnx7oVu/view?usp=sharing, https://drive.google.com/file/d/1dqNY6hz00eb3IycqV98ty--dtXeMv0Zz/view?usp=sharing, https://drive.google.com/file/d/1uB5-htd3rMrX2crVD-j3afYUaDydlwBI/view?usp=sharing, https://drive.google.com/file/d/1vfYmUE9RrTSzsjWVMF6c7Ky30GIuAlXQ/view?usp=sharing, https://drive.google.com/file/d/1vivg83T5Tuol5WAtBzObUX3RwjzE_x1D/view?usp=sharing, https://drive.google.com/file/d/1wrv4uV3zOxkvmsxCutBacI_VEg4dP_eg/view?usp=sharing, https://drive.google.com/file/d/1zYI331I8SuOrQmtxXgR00SpHXIIvF00B/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done

# ## ATAC for monkey Memphis
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/ATAC/Memphis
# LINKFILE='https://drive.google.com/file/d/14PqeuS4sIK65oDouyGX4Kpbf1LfWR2Jp/view?usp=sharing, https://drive.google.com/file/d/169GgOpU47sNttlyZdfu-1L3KPD6_HsHz/view?usp=sharing, https://drive.google.com/file/d/16oKpfpwghYPvitGKWFqzKXW7Vz9bn1yz/view?usp=sharing, https://drive.google.com/file/d/19UxVGZFxYlWRwg4cY5O3Rj4PqreXxghM/view?usp=sharing, https://drive.google.com/file/d/1GdGr9ClvrLdfGrvzTWqNpjuf19M63NRg/view?usp=sharing, https://drive.google.com/file/d/1NVC01YhImqxUbjAHFW8-wPgHHgouuxiI/view?usp=sharing, https://drive.google.com/file/d/1RKZx4Rt-Ld5K9aWpS5jkiUq_VEeWdtky/view?usp=sharing, https://drive.google.com/file/d/1T3I03tNIFUm2x9aouCkO2OQNodxgFOhB/view?usp=sharing, https://drive.google.com/file/d/1TVyMyKRTPcmUh3kR8l_Bbn68aAcIdpQJ/view?usp=sharing, https://drive.google.com/file/d/1Xoqs6yRqCItH9Hi4kdmqaAg3gwBDQIq7/view?usp=sharing, https://drive.google.com/file/d/1fqdTwNzWn6L7mUrffOxfk25ymXQCHiEI/view?usp=sharing, https://drive.google.com/file/d/1fzsxTBnjU3ONU-D3f-Bu4jyruPu0x7o5/view?usp=sharing, https://drive.google.com/file/d/1kUYj_fKe5orE_ssykM3r0X39J7F7SLgP/view?usp=sharing, https://drive.google.com/file/d/1pJu78ZGJtCwNL8Kc3wDofR1j29Bfs3hK/view?usp=sharing, https://drive.google.com/file/d/1tm2_yAsuLhbE6ewXE4ZAvGOpglbass7k/view?usp=sharing, https://drive.google.com/file/d/1z0wPTIMnbApQRAMvaUjDMXJGOQRbtYi1/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done

# ## ATAC for Fiona
# mkdir -p /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/ATAC/Fiona
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/ATAC/Fiona
# LINKFILE='https://drive.google.com/file/d/1RRk3zgllGpWozNI7s8wMJbawMy2136Km/view?usp=sharing, https://drive.google.com/file/d/1qQf8wbdrFJs5c5lIoh7gBWMg7dyeNuC-/view?usp=sharing, https://drive.google.com/file/d/1sninGrUEDxnTzIaI_Se4bDKo7_AsyXmk/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done


# ## RNA for monkeys Oskar and Salem
# LINKFILE='https://drive.google.com/file/d/10GejQt0lQoxFrT8TJbkYqFXdBeo3YWtF/view?usp=sharing, https://drive.google.com/file/d/1IxlImw7iw2aK2tyd0kgVyMjnAP13qXYo/view?usp=sharing, https://drive.google.com/file/d/1jFwFvHtSDo1jEFcS76maoOnngA6W8xLB/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done


# ## RNA for monkey London
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/RNA/London
# LINKFILE='https://drive.google.com/file/d/1-SqvC7peMFmrKDlF5QXGvmh6S3Lx8Hu3/view?usp=sharing, https://drive.google.com/file/d/10WFhb0eIV1N9RY0h4Q332fvTZ8bCc4S9/view?usp=sharing, https://drive.google.com/file/d/14nROADH3bKFT0ROLBcRdolzs4wzrb2lv/view?usp=sharing, https://drive.google.com/file/d/18uiqRdrHfhJe7SAf_i7jqOQq2B3IikOo/view?usp=sharing, https://drive.google.com/file/d/1M2nnK_EC4dszlyKuWx88Jel6swJgdaY0/view?usp=sharing, https://drive.google.com/file/d/1PSGYrry3qqpYfjnlkvBRf7AB4I_Ft6XL/view?usp=sharing, https://drive.google.com/file/d/1QD00jY5WIlqpd3kA19la1P5a7bMyKxR3/view?usp=sharing, https://drive.google.com/file/d/1aIoRJs4b7drOrldfnFPEBZZFxnsR5Ttw/view?usp=sharing, https://drive.google.com/file/d/1hw70x7DzXhDCkfvwtvjKbkgqmlXySvIF/view?usp=sharing, https://drive.google.com/file/d/1idhBhoeQMhv0hcVWOyPvIApgOXqQBu8V/view?usp=sharing, https://drive.google.com/file/d/1ls-1Et69i4qiVWomvWVR74j3brIzUOD3/view?usp=sharing, https://drive.google.com/file/d/1mtIitjmepVG5J5wZJ6Pg1o1nkVBOifRX/view?usp=sharing, https://drive.google.com/file/d/1nDwAgDRXBtXScK_VAXZAeNznq4Xy7AAW/view?usp=sharing, https://drive.google.com/file/d/1oGFEhu6lQb9d5-F2rCtmYQ-pGgYWa6v8/view?usp=sharing, https://drive.google.com/file/d/1oHI9jNSgtrtMyEowHJ6HEkWax0BjRZSG/view?usp=sharing, https://drive.google.com/file/d/1y7RzWzpEJVPbkxfEmTNeI2GAmQSX98zk/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done

# ## RNA for monkey Memphis
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/RNA/Memphis
# LINKFILE='https://drive.google.com/file/d/1-IXtqswrzG3711yqdn0cWtZZMQ0VM3oo/view?usp=sharing, https://drive.google.com/file/d/1-XXeU8qUENKzd6cRsIbthN_Yt6aN1NxT/view?usp=sharing, https://drive.google.com/file/d/11ZcpawMMBcmu6Swchx_li3B06BdJQH8R/view?usp=sharing, https://drive.google.com/file/d/19c2RiwgagDqa-jk7vlN7e-YbAbwlFvjR/view?usp=sharing, https://drive.google.com/file/d/1DpVMM4pkdXhPk2LncxY0WrWqx8XT5AwP/view?usp=sharing, https://drive.google.com/file/d/1Ex53u3S9zJd8q2mVbnBHAXS-PTaZas_h/view?usp=sharing, https://drive.google.com/file/d/1GxCdMO89ek8OjMyAqmcXjBZyPRuWI9Yn/view?usp=sharing, https://drive.google.com/file/d/1JINF8iG4isaF4dimzrQz_K5oHC6kJjcy/view?usp=sharing, https://drive.google.com/file/d/1SIWQDbBivPGuRLi329z1peMnOW8q8KrW/view?usp=sharing, https://drive.google.com/file/d/1Z6wckrqUMI_NYLSPg2PV2wnxQ2p-pYba/view?usp=sharing, https://drive.google.com/file/d/1ZolqH3MzjzhqV8MyVol9WrpUdoITfgYi/view?usp=sharing, https://drive.google.com/file/d/1aXPt4SEtaMsvP9KRuyBkGvf3e49T5aXu/view?usp=sharing, https://drive.google.com/file/d/1c9nzOiJLhozVY_q3oHxACpvwH8SQqpGv/view?usp=sharing, https://drive.google.com/file/d/1dMrL0PSCbQez0py-N5lrJM-E-r7P1ivQ/view?usp=sharing, https://drive.google.com/file/d/1rn2MraL0BXu3os9PhZFKycSX3CWxfS0x/view?usp=sharing, https://drive.google.com/file/d/1tXt1PgNYpD5fcWTwNzFpvxTCs5nQD7l9/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done


# ## RNA for monkey London
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/RNA/London
# LINKFILE='https://drive.google.com/file/d/1PSGYrry3qqpYfjnlkvBRf7AB4I_Ft6XL/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done

# ## RNA for monkey Memphis
# cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/RNA/Memphis
# LINKFILE='https://drive.google.com/file/d/1JINF8iG4isaF4dimzrQz_K5oHC6kJjcy/view?usp=sharing, https://drive.google.com/file/d/1SIWQDbBivPGuRLi329z1peMnOW8q8KrW/view?usp=sharing, https://drive.google.com/file/d/1Z6wckrqUMI_NYLSPg2PV2wnxQ2p-pYba/view?usp=sharing, https://drive.google.com/file/d/1dMrL0PSCbQez0py-N5lrJM-E-r7P1ivQ/view?usp=sharing'
# IDS=$(echo $LINKFILE | tr ',' '\n' | sed 's/^.*\/d\///; s/\/view.*$//')
# for ID in $IDS; do gdown --id ${ID}; done




