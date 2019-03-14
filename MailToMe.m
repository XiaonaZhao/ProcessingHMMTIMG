function MailToMe(inboxAddress)
% https://www.jianshu.com/p/409d2d527326

outboxAddress = 'gittest@163.com';
password = '163gittest';

setpref('Internet','E_mail',outboxAddress);
setpref('Internet','SMTP_Server','smtp.163.com');
setpref('Internet','SMTP_Username',outboxAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','false');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.fallback','false');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(inboxAddress, 'Test passed', 'Maybe it is ok.');
