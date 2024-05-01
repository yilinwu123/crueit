
library(shiny)
library(plotly)
library(DT)
library(ggplot2)
library(tibble)
library(ggsurvfit)
library(survminer)
library(gridExtra)
library(tidyverse)
library(broom)
library(survival)
library(shinydashboard)

select=read.csv("data/subtype_select.csv")
k=4
m=select$cancer[k]
n=select$subtype[k]



data1=read.csv(paste( "data/GDCdata_miRNA-1/",m,"_", n, ".csv",sep=""))
data1<-data1[order(data1$BatchId),]

data1<-data1[,-1]
data=data1[data1$os_time!=0&!is.na(data1$os_time),]
data=data[data$os_time>0,]
names=colnames(data)[1883:ncol(data)]
data2=data[,2:1882]
data3<-log2(data2+1)
genes<-colnames(data3)
datause<-cbind(data3,data$os_status,data$os_time)
colnames(datause)[1882]<-"os_status"
colnames(datause)[1883]<-"os_time"

datausepre=datause[,1:1881]
geneid=colnames(datausepre)
exclude=c()
means1=colMeans(datausepre)
ps_rho_cutoff=0.9
thr=2
x_sub=datausepre[,means1>thr]
xcorr=cor(x_sub)

for(i in 1:ncol(x_sub)){
  for(j in i:ncol(x_sub)){
    if(xcorr[i,j]>ps_rho_cutoff & xcorr[i,j]!=1){
      exclude_1=ifelse(means1[geneid==rownames(xcorr)[i]]>=means1[geneid==colnames(xcorr)[j]],
                       colnames(xcorr)[j], rownames(xcorr)[i])
      exclude=c(exclude, exclude_1)
    }
  }
}
if(length(unique(exclude))==0){
  datapre=x_sub
}else{
  datapre=x_sub[,!(colnames(x_sub) %in% unique(exclude))]
}

geneidfinal=colnames(datapre)
#print(dim(datapre))
#data.train=cbind(data.train.raw$submitter_id,x.train,data.train.raw[,1883:ncol(data.train.raw)])
datausefinal=cbind(datause$os_time,datause$os_status,datapre)
colnames(datausefinal)[1]='os_time'

colnames(datausefinal)[2]='os_status'



#dataplot=read.csv(paste('/users/yilin/Downloads/cureitsurvival-2/cure',m,n,'.csv'))
dataplot=read.csv(paste('/users/yilin/Downloads/cureitsurvival/cure',m,n,'.csv'))
dataplot$color[dataplot$scoretest<0.01]=1
dataplot$color[dataplot$scoretest>=0.01]=0
dataplot$color=as.factor(dataplot$color)
dataplot$curepvalue<-apply(dataplot[,c('curepvalue1','curepvalue2','curepvalue3')],1,min)
dataplot<-dataplot%>%filter(curepvalue!=0)
dataplot$curepvalue=-log10(dataplot$curepvalue)
dataplot$survivalpvalue<-apply(dataplot[,c('survivalpvalue1','survivalpvalue2','survivalpvalue3')],1,min)
dataplot<-dataplot%>%filter(survivalpvalue!=0)
dataplot$survivalpvalue=-log10(dataplot$survivalpvalue)
dataplot$coxpvalue=-log10(dataplot$scoretest)


dataplot<-dataplot[!dataplot$survivalpvalue>=10,]
dataplot<-dataplot[!dataplot$curepvalue>=10,]
dataplot<-dataplot[!dataplot$coxpvalue>=10,]
#dataplot<-dataplot[!abs(dataplot$survival)>=10,]
#dataplot<-dataplot[!abs(dataplot$cure)>=10,]
#dataplot<-dataplot[!abs(dataplot$cox)>=10,]



ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(),
  dashboardBody(
    fluidPage(
      fluidRow(
        column(6,plotlyOutput("plot")),
        column(6,plotlyOutput("text"),plotOutput("text1"))
      )
      #plotlyOutput("plot"),
      #plotlyOutput("text"),
      #plotlyOutput("text1")
    )
  )
)
server<-function(input,output){
  
  #data <- reactive({
  #dataplot
  # })
  
  output$plot <- renderPlotly({
    #d <- data()
    plot_ly(dataplot, x= ~curepvalue, y=~survivalpvalue, split=~color, mode = "markers", type = "scatter", source="mysource")
    #ggplot(data=dataplot,aes(x=curepvalue,y=survivalpvalue,col=color))+geom_point()+scale_color_manual(values=c('Black','Red'))
  })
  
  
  
  
  output$text <- renderPlotly({
    event.data <- event_data("plotly_click", source = "mysource")
    print(event.data)
    if(is.null(event.data)) { return(NULL)}
    index=which(dataplot$survivalpvalue==event.data$y)
    geneid=dataplot$gene[index]
    datausefinal$label1=cut(datausefinal[,geneid],breaks=c(min(datausefinal[,geneid]),quantile(datausefinal[,geneid],0.25),quantile(datausefinal[,geneid],0.5),quantile(datausefinal[,geneid],0.75),max(datausefinal[,geneid])),
                            include.lowest = T,right = F,label=c(0,1,2,3))
    fit<-do.call(survfit,list(Surv(datausefinal$os_time,datausefinal$os_status)~datausefinal$label1,data=datausefinal))
    #plot_ly(datausefinal,x=~os_time,y=~os_status)
    #ggplot(data=datausefinal,aes(x=os_time,y=os_status))+geom_point()
    p<-ggsurvplot(fit,data=datausefinal, risk.table = TRUE, risk.table.y.text.col = TRUE)+labs(title=geneid)
    ggplotly(p[[1]])
  })
  
  
  output$text1 <- renderPlot({
    event.data <- event_data("plotly_click", source = "mysource")
    print(event.data)
    if(is.null(event.data)) { return(NULL)}
    index=which(dataplot$survivalpvalue==event.data$y)
    geneid=dataplot$gene[index]
    #plot_ly(data=datausefinal,x=~str(geneid),type='histogram')
    hist(datausefinal[,geneid],main=geneid)
  })
  
  output$download_plotly_widget <- downloadHandler(
    filename = function() {
      paste("/users/yilin/Downloads/data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # export plotly html widget as a temp file to download.
      saveWidget(as_widget(session_store$plt), file, selfcontained = TRUE)
    }
  )
  
}

shinyApp(ui, server)
