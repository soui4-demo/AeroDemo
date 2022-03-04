#include "StdAfx.h"
#include "SAeroWindow.h"
#include "GaussianBlur.h"

SNSBEGIN

SAeroWindow::SAeroWindow(void):m_radius(5)
{
	m_bLayeredWindow=TRUE;
}

SAeroWindow::~SAeroWindow(void)
{
}

void SAeroWindow::OnCommitSurface(IRenderTarget *pRtDest,LPCRECT pRcDest,IRenderTarget *pRtSrc,LPCRECT pRcSrc,BYTE alpha)
{
	if(!m_midRt)
	{
		__baseCls::OnCommitSurface(pRtDest,pRcDest,pRtSrc,pRcSrc,alpha);
		return;
	}
	CRect rc(pRcDest);
	rc.MoveToXY(0,0);
	m_midRt->BitBlt(&rc,pRtDest,pRcDest->left,pRcDest->top,SRCCOPY);
	m_midRt->AlphaBlend(&rc,pRtSrc,pRcSrc,alpha);
	IBitmapS * pBmp = (IBitmapS*)m_midRt->GetCurrentObject(OT_BITMAP);
	BYTE *pData = (BYTE*)pBmp->LockPixelBits();
	SASSERT(LINEWID(rc.Width()*4)==rc.Width()*4);
	memcpy(m_buf.ptr(),pData,rc.Height()*rc.Width()*4);
	CGaussianBlur::GaussianBlur4(m_buf.ptr(),pData,rc.Width(),rc.Height(),4,m_radius);
	pBmp->UnlockPixelBits(pData);
	pRtDest->BitBlt(pRcDest,m_midRt,0,0,SRCCOPY);
}

void SAeroWindow::OnSize(UINT nType, CSize size)
{
	__baseCls::OnSize(nType,size);
	CRect rcWnd = GetWindowRect();
	m_midRt = NULL;
	GETRENDERFACTORY->CreateRenderTarget(&m_midRt,rcWnd.Width(),rcWnd.Height());
	m_buf = new BYTE[rcWnd.Height()*LINEWID(rcWnd.Width()*4)];
}

SNSEND